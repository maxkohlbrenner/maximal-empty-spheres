import numpy as np
from lie_geometry.lie_geometry import lie_inner_product, circle_to_lie, lie_to_circle, plane_to_lie
from lie_geometry.lie_geometry import line_quadric_intersection as line_quadric_intersection, line_quadric_intersection_single
from scipy.spatial import ConvexHull
from numba import jit
import time

from gpytoolbox import point_cloud_to_mesh

class LieConeSDFReconstruction:

    def __init__(self, G, N_extra=None, filter_type=0, n_check_planes=10, cut_bbx_factor=-1., timeit=False, verbose=True,filter_results=False,psr_screening_weight=1.):

        self.timings = []
        if filter_type == 3:
            # combine pos+neg points
            filter_types = [1,2]
            cut_bbx_factors = [cut_bbx_factor,-1]
        else:
            filter_types = [filter_type]
            cut_bbx_factors = [cut_bbx_factor]

        if verbose:
            print("Run cone construction(s)")
        self.lcs   = [LieCone(G, N_extra, ft,
        n_check_planes=n_check_planes,
        cut_bbx_factor=cbf,
        timeit=timeit,
        filter_results=filter_results) for ft,cbf in zip(filter_types,cut_bbx_factors)]

        if verbose:
            print("Extract closest points")

        self.psets = [lc.calc_psr_contact_points() for lc in self.lcs]
        self.points  = np.concatenate([P for (P,N) in self.psets],axis=0)
        self.normals = np.concatenate([N for (P,N) in self.psets],axis=0)

        for lc in self.lcs:
            if len(lc.timings) > 0:
                for (tc,t) in lc.timings:
                    print(f"{tc:<30}: {t:1.4}s")


        if verbose:
            print("Run PSR")
        self.V, self.F = point_cloud_to_mesh(self.points, self.normals, psr_screening_weight=psr_screening_weight)

class LieCone:

    def __init__(self, G, N_extra=None, filter_type=0, n_check_planes=10, cut_bbx_factor=-1., timeit=False, verbose=False, rot_around_xdp2=False, filter_results=False, assertions=False, full_simplices_only=True, atol=1e-8):

        self.timings = []
        self.gather_timings = timeit
        if timeit:
            strt_method = time.time()
            strt = strt_method

        self.G = G.copy()
        self.D = self.G.shape[1]-1
        self.n_check_planes=n_check_planes

        self.H = lie_inner_product(self.D)

        # filter out certain spheres
        if verbose:
            print(" ### Preprocessing ###")
        if filter_type == 1:
            # POS
            G_filtered = G[G[:,-1] >= 0.]
            self.signs = np.ones(G_filtered.shape[0])
        elif filter_type == 2:
            # NEG
            G_filtered = G[G[:,-1]  < 0.]
            self.signs = - np.ones(G_filtered.shape[0])
        else:
            # BOTH
            G_filtered = G
            self.signs = np.sign(G_filtered[:,-1])

        # pos only
        SE = G_filtered.copy()
        SE[:,-1] = np.abs(SE[:,-1])
        self.SE = SE

        if assertions:
            assert np.all(self.SE[:,-1] >= atol), "All input radii must be positive"
        if timeit:
            stp = time.time()
            self.timings.append(("Preprocessing", stp-strt))
            strt = stp

        if verbose:
            print("... transform to lie plane normals")

        SL = circle_to_lie(self.SE)
        NC = SL@self.H

        if assertions:
            assert np.allclose(np.einsum("ni,ij,nj->n",SL,self.H,SL),0.), "All spheres must lie on the Lie quadric"

        if N_extra is not None:
            NC = np.concatenate([S,N_extra],axis=0)

        # define rotation
        s0 = np.zeros(self.D+3)
        s0[-1] = -1
        s0[-3] = -1
        dst = np.zeros(self.D+3)
        dst[-1] = 1.
        R = self.rotation_from_to(s0,dst)

        NCR = NC@R        # Rotation
        NCR /= NCR[:,-1:] # Normalize
        Np = NCR[:,:-1]   # Intersect

        if timeit:
            stp = time.time()
            self.timings.append(("Normal Creation", stp-strt))
            strt = stp

        # Convex Hull Computation (MAIN)
        if verbose:
            print("... calculate convex hull")
        #self.ch = ConvexHull(N)
        self.ch = ConvexHull(Np, qhull_options='Q12')

        if timeit:
            stp = time.time()
            self.timings.append(("Convex Hull", stp-strt))
            strt = stp

        #calculate kernel per simplex
        if verbose:
            print("... calculate kernel per simplex")
        SC = NC[self.ch.simplices]
        U_, s_, Vt_ = np.linalg.svd(SC)
        Ks = Vt_[:,-1]

        K_is_full_simplex_mask = np.min(np.absolute(s_),axis=1) > atol

        if verbose:
            dots = np.einsum('njk,nk->nj', SC,Ks)
            print("    ->Kernels valid (per simplex): {}".format(np.allclose(dots, 1e-15)))

        Kds_ = Ks@NC[:self.n_check_planes].T
        inv = Kds_.max(axis=1) >= 1e-5
        Ks[inv] *= -1
        self.Ks = Ks

        if assertions:
            if full_simplices_only:
                print("☠️ Only consider K simplices with full rank. ☠️")
                assert np.all(Ks[K_is_full_simplex_mask]@NC.T<=atol), "All directions d_i must lie in the cone"
            else:
                print("☠️☠️☠️ DEACTIVATED K CHECK, WERE ALSO USING NON-FULL SIMPLICES ☠️☠️☠️")
                assert np.all(Ks@NC.T<=atol), "All directions d_i must lie in the cone"

        if timeit:
            stp = time.time()
            self.timings.append(("HS Intersection", stp-strt))
            strt = stp

        if verbose:
            print("... calculate edge and facet connectivity")
        self.Es, self.E_ninds, self.Ts = self.construct_edges_and_facets(self.ch, True)

        if timeit:
            stp = time.time()
            self.timings.append(("Connectivity", stp-strt))
            strt = stp

        if verbose:
            print("... find solutions on the edges")

        if full_simplices_only:
            print("☠️ Only consider edges connecting two full simplices ☠️")
            Es_between_full_simplices = np.all(K_is_full_simplex_mask[self.Es],axis=1)
            self.solutions, self.T_indices = line_quadric_intersection(self.Ks[self.Es[Es_between_full_simplices]],self.H)
        else:
            self.solutions, self.T_indices = line_quadric_intersection(self.Ks[self.Es],self.H)

        if timeit:
            stp = time.time()
            self.timings.append(("Lie Intersection", stp-strt))
            self.timings.append(("Total", stp-strt_method))
            strt = stp

        if verbose:
            print("... DONE\n")

        if timeit:
            for (tc,t) in self.timings:
                print(f"{tc:<30}: {t:1.4}s")

        if cut_bbx_factor > 0.:
            # cut planes
            P = self.bbx_planes(cut_bbx_factor)
            #valid_solutions = ((solutions@S.T).max(axis=1) <= 1e-10)
            print(self.solutions.shape, self.H.shape, P.shape)
            cp_flags = ((self.solutions@self.H@P.T) < 0.).all(axis=1)
        else:
            cp_flags = True

        # finite radius flag: solution sometimes contains very spheres becuse radius comp. went near zero
        fr = np.abs(self.solutions[:,-1]) >= 1e-10

        # negative convention and negative radius flags
        nc = self.solutions[:,-1] < 0.
        nr = np.sign(self.solutions[:,-1]) != np.sign(self.solutions[:,-2])

        if filter_results:
            # highly illegal: just post_filter result
            print("-----------> ☠️☠️☠️☠️POST FILTER ☠️☠️☠️ <-----------")
            fint = 1000
            If = np.random.choice(np.arange(NC.shape[0]), fint)
            resultfilter = (self.solutions@NC[If].T <= 1e-10).all(axis=1)
        else:
            resultfilter = 1

        self.keep_mask = (nc*nr*fr*cp_flags*resultfilter).astype("bool")
        self.selected_solutions = self.solutions[self.keep_mask]

        self.sol = lie_to_circle(self.selected_solutions)
        self.sol_T_indices = self.T_indices[self.keep_mask]
        if full_simplices_only:
            self.sol_Ts = self.Ts[Es_between_full_simplices][self.sol_T_indices]
        else:
            self.sol_Ts = self.Ts[self.sol_T_indices]

        #T_centers = (SE[:,:-1][Ts]).mean(axis=1)
        #sol_T_centers = T_centers[sol_T_indices]



    def calc_psr_contact_points(self):

        if self.gather_timings:
            strt = time.time()
        self.ami  = LieCone.get_maximal_adjacent_spheres(self.SE, self.sol, self.sol_Ts)
        if self.gather_timings:
            stp = time.time()
            self.timings.append(("Maximal Adjacent Spheres", stp-strt))
            strt = stp

        Pts = self.get_psr_oriented_points()

        if self.gather_timings:
            stp = time.time()
            self.timings.append(("Psr Points", stp-strt))
            strt = stp
        return Pts

    def bbx_planes(self, cf):
        return LieCone.bbx_planes_(self.SE,cf)

    @staticmethod
    def bbx_planes_(G, factor):
        Dim = G.shape[1]-1
        A = np.concatenate([np.eye(Dim),-np.eye(Dim)],axis=0)
        b = factor * np.concatenate([-np.min(G[:,:-1],axis=0), np.max(G[:,:-1],axis=0)],axis=0)
        P = plane_to_lie(np.concatenate([A,b[:,None]],axis=1))
        return P



    @staticmethod
    @jit(nopython=True)
    def get_maximal_adjacent_spheres(SE, sol, sol_Ts):
        sol_rads = sol[:,-1]
        adjacent_minrads = np.max(np.absolute(sol_rads)*1.1) * np.ones(SE.shape[0])
        adjacent_minrad_indices = - np.ones(SE.shape[0],dtype="int")

        for si, (sr, sTi) in enumerate(zip(sol_rads,sol_Ts)):
            repl_inds = sTi[(adjacent_minrads[sTi] > sr)]
            adjacent_minrads[repl_inds] = sr
            adjacent_minrad_indices[repl_inds] = si

        return adjacent_minrad_indices

    #@jit
    def get_psr_oriented_points(self):

        #self.ami = self.get_maximal_adjacent_spheres()

        ami_valid = self.ami >= 0

        N =  self.SE[ami_valid,:-1] - self.sol[self.ami[ami_valid]][:,:-1]
        N /= np.linalg.norm(N,axis=1,keepdims=True)
        P = self.SE[ami_valid,:-1] - np.abs(self.SE[ami_valid,-1:])*N
        N *= self.signs[ami_valid,None]

        self.P = P
        self.N = N

        return P, N

    def solutions_on_edges(self, Es, Ks):
        solutions = []
        T_indices = []

        for e_idx,e in enumerate(Es):
            k = Ks[e[0]]
            n = Ks[e[1]]
            for l in line_quadric_intersection_single(k,n, self.H):
                    if 0 <= l <= 1.:
                        solutions.append((1-l)*k+l*n)
                        T_indices.append(e_idx)

        solutions = np.array(solutions)
        T_indices = np.array(T_indices)

        if False:
            print("Found {} solution".format(solutions.shape[0]))
            print("Solutions on the lie quadric: {}".format(np.allclose(np.einsum("ni,ij,nj->n",solutions,H,solutions), 0.)))
            #print("Maximum over the edge value: {}".format((solutions@S.T).max(axis=1)))

        return solutions, T_indices

    def construct_edges_and_facets(self, ch, add_facet_vertices=False):
        F = ch.neighbors
        Fi = np.arange(F.shape[0])
        Fi_b = np.broadcast_to(Fi[:,None], F.shape)
        #Hes = np.stack([F, Fi_b],axis=-1).reshape(-1,2)
        Hes = np.stack([Fi_b, F],axis=-1).reshape(-1,2)
        pos_edges = Hes[:,0] < Hes[:,1]

        # Es: nedges, 2
        # gives the two simplex indices of an edge, in ascending order
        Es = Hes[pos_edges]

        # E_ninds: nedges, 2
        # for each edge in E: gives the first simplex index and the index of the face in the simplces
        # (use to extract the common vertices)
        Ni = np.arange(F.shape[1])
        Ni_b = np.broadcast_to(Ni[None,:], F.shape)
        He_ninds = np.stack([Fi_b, Ni_b],axis=-1).reshape(-1,2)
        E_ninds = He_ninds[pos_edges]

        # check: reconstruct edges from Edge neighbor indices:
        if False:
            E__ = np.stack([E_ninds[:,0], F[E_ninds[:,0],E_ninds[:,1]]],axis=1)
            print("Edge neighbor indices consistent? : {}".format(np.allclose(Es, E__)))

        if add_facet_vertices:

            # an edge connects two simplices along a face
            # Ts are the vertices of this face
            #Ts = (F[E_ninds[:,0]])[(E_ninds[:,1,None] != Ni[None,:])].reshape(-1,Ni_b.shape[1]-1)
            facet_vertex_indices = (E_ninds[:,1,None] != Ni[None,:])
            Ts = (ch.simplices[E_ninds[:,0]])[facet_vertex_indices].reshape(-1,Ni_b.shape[1]-1)

            if False:
                print("Facets all contained in both Edge simplices? {}".format((Ts[:,None,None,:] == ch.simplices[Es][:,:,:,None]).any(axis=(2,3)).all(axis=1).all()))

            return Es, E_ninds, Ts

        else:
            return Es, E_ninds

    def load_results(self):
        pass

    @staticmethod
    def rotation_from_to(source, target):
        D = source.shape[0]
        S = np.zeros((D,D))

        S[:,0] = source
        U,s,Vt = np.linalg.svd(S, full_matrices=True)

        Rs = U@Vt
        if np.linalg.det(Rs) < 0.:
            Rs *= -1.

        T = np.zeros((D,D))
        T[:,0] = target

        U,S,Vt = np.linalg.svd(T)
        Rt = U@Vt
        if np.linalg.det(Rt) < 0.:
            Rt *= -1.

        return Rt@Rs.T
