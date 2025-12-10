import os
from os.path import join
import subprocess
import numpy as np
from gpytoolbox import point_cloud_to_mesh

basefolder = join(os.path.dirname(os.path.abspath(__file__)))

def MESReconstruction(G,D,screening_weight=1., return_oriented_points=False, cleanup=True, save_folder="./"):

    cpp_bin_path = join(basefolder, "build/empty_spheres_reconstruction")
    if not os.path.exists(cpp_bin_path):
        print("ERROR (maximal empty spheres reconstruction:) CPP executable NOT FOUND!")
        return False
    else:

        gridpath = join(save_folder,"tmp.csv")
        np.savetxt(gridpath, np.concatenate([G,D[:,None]],axis=1), delimiter=",")

        subprocess.run([cpp_bin_path, gridpath])
        pwnpath = './pwn.csv'
        cgal_pwns = np.genfromtxt(pwnpath, delimiter=',')

        if cleanup:
            for pth in [gridpath,pwnpath]:
                os.remove(pth)

        R = point_cloud_to_mesh(cgal_pwns[:,:3],cgal_pwns[:,3:],psr_screening_weight=screening_weight)

        if return_oriented_points:
            return *R, cgal_pwns[:,:3], cgal_pwns[:,3:]
        else:
            return R
