import numpy as np
import sys
sys.path.append('/Users/denis/ownCloud/pfc-layers/analysis/fmri-analysis/library')
sys.path.append('/data/p_02389/code/fmri-analysis/library')
import voxeldepths_from_surfaces as vfs
import nibabel as nib
import os
from numba import jit, prange
from skimage.segmentation import expand_labels
import matplotlib.pyplot as plt
from joblib import Parallel, delayed, parallel_backend

@jit(nopython=True)
def surf_inside_out_f(M,V,xyz):
    """Surf inside out single point version, to be used in conjunction with joblib.Parallel.
    """
    
    n_tris = M.shape[2]
    n = 0
    p_x, p_y, p_z = xyz

    for i_tris in range(n_tris):
        # assume pc_x = 1, pc_y = 0, pc_z = 0
        pc_n, pc_s, pc_t = M[0,:,i_tris]

        vp_x = p_x - V[i_tris,0]
        vp_y = p_y - V[i_tris,1]
        vp_z = p_z - V[i_tris,2]

        p_n = vp_x*M[0,0,i_tris] + vp_y*M[1,0,i_tris] + vp_z*M[2,0,i_tris]
        p_s = vp_x*M[0,1,i_tris] + vp_y*M[1,1,i_tris] + vp_z*M[2,1,i_tris]
        p_t = vp_x*M[0,2,i_tris] + vp_y*M[1,2,i_tris] + vp_z*M[2,2,i_tris]
            
        alpha = p_n * pc_s
        beta =  p_s * pc_n
        gamma = p_n * pc_t
        delta = p_t * pc_n

        if p_n >= 0:
            if alpha>=beta and gamma>=delta and alpha+gamma<=beta+delta-pc_n:
                n = n + 1
        elif alpha<=beta and gamma<=delta and alpha+gamma>=beta+delta-pc_n:
            n = n + 1

    return n

@jit(nopython=True)
def surf_inside_out_loop(M,V,XYZ):
    """Reimplementation of mex function  dcSurfInsideOutLoop(M,V,XYZ).c
    part of pipeline to determine gm,wm,outside for arbitrary point based on surfaces,
    and subsequent cortical depth calculation
    M - list of 3x3 matriced
    V - list of 3d vectors, base vertices
    XYZ - list of 3d vectors, points to test
    """
    
    n_tris = M.shape[2]
    n_voxels = XYZ.shape[0]
    n = np.zeros(n_voxels)

    # loop over triangles
    for i_tris in range(n_tris):
        v_x, v_y, v_z = V[i_tris,:]
        M0, M1, M2 = M[0,:,i_tris]
        M3, M4, M5 = M[1,:,i_tris]
        M6, M7, M8 = M[2,:,i_tris]

        # assume pc_x = 1, pc_y = 0, pc_z = 0
        pc_n = M0
        pc_s = M1
        pc_t = M2

        # loop over voxels
        for i_voxel in range(n_voxels):
            p_x, p_y, p_z = XYZ[i_voxel,:]

            vp_x = p_x - v_x
            vp_y = p_y - v_y
            vp_z = p_z - v_z

            p_n = vp_x*M0 + vp_y*M3 + vp_z*M6
            p_s = vp_x*M1 + vp_y*M4 + vp_z*M7
            p_t = vp_x*M2 + vp_y*M5 + vp_z*M8
            
            alpha = p_n * pc_s
            beta =  p_s * pc_n
            gamma = p_n * pc_t
            delta = p_t * pc_n

            if p_n >= 0:
                if alpha>=beta and gamma>=delta and alpha+gamma<=beta+delta-pc_n:
                    n[i_voxel] = n[i_voxel] + 1
            elif alpha<=beta and gamma<=delta and alpha+gamma>=beta+delta-pc_n:
                n[i_voxel] = n[i_voxel] + 1
    return n
        
@jit(nopython=True)
def surf_bounding_voxels(A,B,C,array_size):
    n_tris = A.shape[0]
    mask = np.zeros(array_size,dtype=np.int0)
    for i_tris in range(n_tris):
        tri_min_x = max(0, np.floor(min(A[i_tris,0],B[i_tris,0],C[i_tris,0])))
        tri_min_y = max(0, np.floor(min(A[i_tris,1],B[i_tris,1],C[i_tris,1])))
        tri_min_z = max(0, np.floor(min(A[i_tris,2],B[i_tris,2],C[i_tris,2])))
        tri_max_x = min(array_size[0]-1,
                        np.ceil(max(A[i_tris,0],B[i_tris,0],C[i_tris,0])))
        tri_max_y = min(array_size[1]-1,
                        np.ceil(max(A[i_tris,1],B[i_tris,1],C[i_tris,1])))
        tri_max_z = min(array_size[2]-1,
                        np.ceil(max(A[i_tris,2],B[i_tris,2],C[i_tris,2])))
        for i_x in range(tri_min_x, tri_max_x+1):
            for i_y in range(tri_min_y, tri_max_y+1):
                for i_z in range(tri_min_z, tri_max_z+1):
                    mask[i_x,i_y,i_z] = 1
    return mask

@jit(nopython=True)
def calc_M(A,B,C):
    n_tris = A.shape[0]
    # triangle edges:
    AB = B - A
    CA = A - C
    # triangle normals
    N = np.cross(AB,-CA)
    # define M
    M = np.zeros((3,3,n_tris))
    for i_tris in range(n_tris):
        M[:,:,i_tris] = np.linalg.inv(np.vstack((N[i_tris,:].T,
                                                 AB[i_tris,:].T,
                                                 -CA[i_tris,:].T)))
    return M

def segment_surf_sides(surf,array_size,gm_ribbon=None,n_jobs=4):
    A = surf['vertices'][surf['tris'][:,0],:]
    B = surf['vertices'][surf['tris'][:,1],:]
    C = surf['vertices'][surf['tris'][:,2],:]
    # calculate GM and WM enclosing voxels
    
    surf_boundary = surf_bounding_voxels(A,B,C,array_size)
    if gm_ribbon is not None:
        surf_boundary = surf_boundary & np.isnan(gm_ribbon)
        
    #  prepare list of voxel coordinates XYZ
    XYZ = np.array(np.argwhere(surf_boundary),dtype=np.float)
    XYZidcs = np.where(surf_boundary)

    M = calc_M(A,B,C)

    with parallel_backend('loky', inner_max_num_threads=1):
        n = Parallel(n_jobs=n_jobs)(
            delayed(
                lambda xyz:
                surf_inside_out_f(M,A,xyz))(xyz)
            for xyz in XYZ)
           
    S = np.zeros(array_size)
    S[XYZidcs] = np.mod(n,2)+1

    return S, surf_boundary
    
def calc_segmentation_from_surf(surf_pial_file,surf_white_file,volume_file,upsample_factor=None,gm_ribbon=None,n_jobs=4):
    volume = nib.load(volume_file)

    n_x, n_y, n_z = volume.shape
    voxel_to_scanner = volume.affine

    if upsample_factor is not None:
        n_x = int(np.floor(n_x * upsample_factor))
        n_y = int(np.floor(n_y * upsample_factor))
        n_z = int(np.floor(n_z * upsample_factor))
        upsampled_to_voxel = np.matrix([[1/upsample_factor,0,0,-1/(2*upsample_factor)],
                                       [0,1/upsample_factor,0,-1/(2*upsample_factor)],
                                       [0,0,1/upsample_factor,-1/(2*upsample_factor)],
                                       [0,0,0,1]]);
        grid_to_scanner = voxel_to_scanner @ upsampled_to_voxel
    else:
        grid_to_scanner = voxel_to_scanner
    
    surf_pial = vfs.load_fs_surf_in_grid(surf_pial_file,grid_to_scanner)
    surf_white = vfs.load_fs_surf_in_grid(surf_white_file,grid_to_scanner)

    seg_white, boundary_white  = segment_surf_sides(surf_white,(n_x,n_y,n_z),
                                                      gm_ribbon=gm_ribbon,n_jobs=n_jobs)
    seg_pial, boundary_pial = segment_surf_sides(surf_pial,(n_x,n_y,n_z),
                                                   gm_ribbon=gm_ribbon,n_jobs=n_jobs)

    # add test for overlaps?
    if gm_ribbon is not None:
        ribbon_boundary_seg = 3 * ~np.isnan(gm_ribbon) + \
            1 * (seg_pial==1) + 2 * (seg_white==2)
    else:
        ribbon_boundary_seg = 3 * (((seg_white==1) | (seg_pial==2)) & ~((seg_pial==1)|(seg_white==2))) + \
            1 * (seg_pial==1) + 2 * (seg_white==2)
    
    ribbon_boundary_seg_expand = expand_labels(ribbon_boundary_seg,max(n_x,n_y,n_z)*2)

    return ribbon_boundary_seg, ribbon_boundary_seg_expand, seg_white, seg_pial, surf_white, surf_pial

def calc_boundary_seg_hemi(surf_pial,surf_white,grid_size,gm_ribbon=None,n_jobs=4):

    boundary_seg_white, boundary_white  = \
        segment_surf_sides(surf_white,grid_size,gm_ribbon=gm_ribbon,n_jobs=n_jobs)
    boundary_seg_pial, boundary_pial  = \
        segment_surf_sides(surf_pial,grid_size,gm_ribbon=gm_ribbon,n_jobs=n_jobs)

    if gm_ribbon is not None:
        boundary_seg_ribbon = 3 * ~np.isnan(gm_ribbon) + \
            1 * (boundary_seg_pial==1) + 2 * (boundary_seg_white==2)
    else:
        boundary_seg_ribbon= 3 * (((boundary_seg_white==1) | (boundary_seg_pial==2)) &
                                  ~((boundary_seg_pial==1)|(boundary_seg_white==2))) + \
                                  1 * (boundary_seg_pial==1) + 2 * (boundary_seg_white==2)
    return boundary_seg_ribbon,boundary_seg_white, boundary_seg_pial,boundary_white, boundary_pial


def calc_segmentation_from_both_hemi_surfs(lh_surf_pial_file,lh_surf_white_file,
                                           rh_surf_pial_file,rh_surf_white_file,
                                           seg_ribbon_fname,
                                           volume_file,upsample_factor=None,gm_ribbon=None,n_jobs=4):
    volume = nib.load(volume_file)
    if type(gm_ribbon) is nib.nifti1.Nifti1Image:
        gm_ribbon = gm_ribbon.get_fdata()
    
    n_x, n_y, n_z = volume.shape
    voxel_to_scanner = volume.affine

    if upsample_factor is not None:
        n_x = int(np.floor(n_x * upsample_factor))
        n_y = int(np.floor(n_y * upsample_factor))
        n_z = int(np.floor(n_z * upsample_factor))
        upsampled_to_voxel = np.matrix([[1/upsample_factor,0,0,-1/(2*upsample_factor)],
                                       [0,1/upsample_factor,0,-1/(2*upsample_factor)],
                                       [0,0,1/upsample_factor,-1/(2*upsample_factor)],
                                       [0,0,0,1]]);
        grid_to_scanner = voxel_to_scanner @ upsampled_to_voxel
    else:
        grid_to_scanner = voxel_to_scanner

    lh_surf_pial = vfs.load_fs_surf_in_grid(lh_surf_pial_file,grid_to_scanner)
    lh_surf_white = vfs.load_fs_surf_in_grid(lh_surf_white_file,grid_to_scanner)
    rh_surf_pial = vfs.load_fs_surf_in_grid(rh_surf_pial_file,grid_to_scanner)
    rh_surf_white = vfs.load_fs_surf_in_grid(rh_surf_white_file,grid_to_scanner)

    xform = grid_to_scanner
    
    lh_boundary_seg_ribbon,lh_boundary_seg_white,lh_boundary_seg_pial,lh_boundary_white, lh_boundary_pial = \
        calc_boundary_seg_hemi(lh_surf_pial, lh_surf_white,(n_x,n_y,n_z),gm_ribbon,n_jobs)
    
    rh_boundary_seg_ribbon,rh_boundary_seg_white,rh_boundary_seg_pial,rh_boundary_white, rh_boundary_pial = \
        calc_boundary_seg_hemi(rh_surf_pial, rh_surf_white,(n_x,n_y,n_z),gm_ribbon,n_jobs)

    boundary_seg_ribbon = 3 * ((lh_boundary_seg_ribbon==3) | (rh_boundary_seg_ribbon==3)) + \
        2 * ((lh_boundary_seg_ribbon==2)| (rh_boundary_seg_ribbon==2)) + \
        1 * (((lh_boundary_seg_ribbon==1) & ~((rh_boundary_seg_ribbon==2)|(rh_boundary_seg_ribbon==3))) | \
             ((rh_boundary_seg_ribbon==1) & ~((lh_boundary_seg_ribbon==2)|(lh_boundary_seg_ribbon==3))))
    nib.save(nib.nifti1.Nifti1Image(boundary_seg_ribbon, xform),
             'boundary_seg_ribbon.nii')
    
    seg_ribbon = expand_labels(boundary_seg_ribbon,max(n_x,n_y,n_z)*2)
    
    xform = grid_to_scanner
    nii_seg_ribbon = nib.nifti1.Nifti1Image(seg_ribbon, xform)
    nib.save(nii_seg_ribbon, seg_ribbon_fname)
    
    return nii_seg_ribbon

def edot(A,B):
    return A[0,:]*B[0,:]+A[2,:]*B[2,:]+A[2,:]*B[2,:]

def plot_surf_zslice(surf,z,c='black',lw=1,transpose=False):
    # define a point and normal vector for slicing plane
    p0 = np.array([0, 0, z])
    n = np.array([0, 0, 1])
    
    # q_[A,B,C] = corner vertices
    # p_[A,B,C] = vectors from one corner to next
    qA = surf['vertices'][surf['tris'][:,0],:].T
    pA = surf['vertices'][surf['tris'][:,1],:].T - qA

    qB = surf['vertices'][surf['tris'][:,1],:].T
    pB = surf['vertices'][surf['tris'][:,2],:].T - qB

    qC = surf['vertices'][surf['tris'][:,2],:].T
    pC = surf['vertices'][surf['tris'][:,0],:].T - qC

    # if triangle edge lines given as q_[A,B,C] + t_[A,B,C] * p_[A,B,C],
    # compute parameter t_[A,B,C] of intersection point with slice plane
    # (t_[A,B,C] = \frac{(p0 - q_[A,B,C]) \cdot n)}{p_[A,B,C] \cdot n}
    tA = edot(p0[:,None]-qA,n[:,None])/edot(pA,n[:,None])
    tB = edot(p0[:,None]-qB,n[:,None])/edot(pB,n[:,None])
    tC = edot(p0[:,None]-qC,n[:,None])/edot(pC,n[:,None])

    # check which edge intersections are between triangle corners
    cutTrisA = (tA>=0) & (tA<=1)
    cutTrisB = (tB>=0) & (tB<=1)
    cutTrisC = (tC>=0) & (tC<=1)
    # check for pairs of edges that intersect with slice
    cutTrisAB = (cutTrisA) & (cutTrisB)
    cutTrisBC = (cutTrisB) & (cutTrisC)
    cutTrisCA = (cutTrisC) & (cutTrisA)

    # calculate intersection points
    aA1 = qA[:,cutTrisAB]+tA[cutTrisAB]*pA[:,cutTrisAB]
    aB1 = qB[:,cutTrisAB]+tB[cutTrisAB]*pB[:,cutTrisAB]
    aB2 = qB[:,cutTrisBC]+tB[cutTrisBC]*pB[:,cutTrisBC]
    aC2 = qC[:,cutTrisBC]+tC[cutTrisBC]*pC[:,cutTrisBC]
    aC3 = qC[:,cutTrisCA]+tC[cutTrisCA]*pC[:,cutTrisCA]
    aA3 = qA[:,cutTrisCA]+tA[cutTrisCA]*pA[:,cutTrisCA]
    
    if transpose == True:
        dims = [1,0]
    else:
        dims = [0,1] 
    for i in range(aA1.shape[1]):
        plt.plot([aA1[dims[0],i], aB1[dims[0],i]],[aA1[dims[1],i], aB1[dims[1],i]],color=c,linewidth=lw)
    for i in range(aB2.shape[1]):
        plt.plot([aB2[dims[0],i], aC2[dims[0],i]],[aB2[dims[1],i], aC2[dims[1],i]],color=c,linewidth=lw)
    for i in range(aC3.shape[1]):
        plt.plot([aC3[dims[0],i], aA3[dims[0],i]],[aC3[dims[1],i], aA3[dims[1],i]],color=c,linewidth=lw)


def plot_one(img,axis,title_str,surfs):
    fig=plt.figure(figsize=(12,8), dpi= 300, facecolor='w', edgecolor='k')
    z = int(np.floor(img.shape[2]/2))
    plt.imshow(img[:,:,z],origin='lower')
    plt.axis(axis)
    plt.title(title_str)
    for surf,c in surfs:
        plot_surf_zslice(surf,z,c=c,lw=2,transpose=True)

def plot_two(img_white,img_pial,axis,title_str,surfs):
    plot_one(img_white,axis,'White matter surface, '+title_str,[surfs[0]])
    plt.show()
    plot_one(img_pial,axis,'Pial surface, '+title_str,[surfs[1]])
    plt.show()   

def visualize(ribbon, intermediate_results,f=1):
   
    axis = [f*40,f*140,f*130,f*160]
    surf_pial, surf_white, intermediate_pial, intermediate_white, S_expand_pial, S_expand_white = intermediate_results
    surf_boundary_pial, S_pial = intermediate_pial
    surf_boundary_white, S_white = intermediate_white
    surfs = ((surf_white,'green'),(surf_pial,'red'))
    
    plot_two(surf_boundary_white*2,surf_boundary_pial*2,axis,'surface bounding voxels',surfs)
    plot_two(S_white,S_pial,axis,'segmentation of bounding voxels',surfs)
    plot_two(S_expand_white,S_expand_pial,axis,
             'morphological expansion of boundary segmentation to entire volume',surfs)
    plot_one(ribbon-1,axis,'Combining both surface segmentation into cortical ribbon segmentation',surfs)
    plt.show()
