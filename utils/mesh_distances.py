import numpy as np
import gpytoolbox as gpy

def mdists(V,F, V_gt, F_gt):
    # sdf values and contact information
    sds, inds, lds = gpy.signed_distance(V, V_gt, F_gt)

    # closest points and sdf normals:
    cps = (lds[:,None,:] @  V_gt[F_gt][inds]).squeeze()

    M = gpy.massmatrix(V, F, type="barycentric")
    ds = np.linalg.norm(V - cps,axis=1)

    # weighted distance errors
    li = ds.max()
    l1 = (M@ds).mean() / M.sum()
    l2 = ds@M@ds / M.sum()

    return li,l1,l2

def ts_distances(R,GT,upsample_iters=1,upsample_gt=False):
    if upsample_iters>0:
        R_  = gpy.subdivide(*R,method='upsample',iters=upsample_iters)[:2]
        
        if upsample_gt:
            GT_ = gpy.subdivide(*GT,method='upsample',iters=upsample_iters)[:2]
        else:
            # only upsample reconstructions (GT assumed to be high res)
            GT_ = GT
    else:
        R_  = R
        GT_ = GT
        
    dto = mdists(*R_, *GT_)
    dfr = mdists(*GT_,*R_) 
    
    # Hausdorff, Mean abs, Chamfer distances
    return np.max([dto[0],dfr[0]]), np.mean([dto[1],dfr[1]]), np.sqrt(dto[2]) + np.sqrt(dfr[2])
