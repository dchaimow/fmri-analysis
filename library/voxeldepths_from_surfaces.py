import numpy as np
import nibabel as nib
from scipy.optimize import root_scalar, minimize_scalar
from numba import jit
from numba.typed import List
from joblib import Parallel, delayed, parallel_backend

def load_fs_surf_in_scanner_space(surf_file):
    surf = load_fs_surf_in_grid(surf_file,np.eye(4))
    return surf

def load_fs_surf_in_grid(surf_file,grid_to_scanner):
    coords, faces, volume_info = nib.freesurfer.read_geometry(surf_file,read_metadata=True)
    surf = dict()
    surf['nVertices'] = coords.shape[0]
    surf['nTris'] = faces.shape[0]
    fs_to_scanner=np.eye(4)
    fs_to_scanner[:3,3]=volume_info['cras']
    grid_to_scanner = np.matrix(grid_to_scanner)
    fs_to_voxel = grid_to_scanner.I @ fs_to_scanner
    surf['vertices'] = np.array(fs_to_voxel.dot(
        np.hstack((coords,np.ones((surf['nVertices'],1)))).T).T[:,:3])
    surf['tris'] = faces
    return surf

def merge_surfaces(surf1, surf2):
    surf=dict()
    surf['nTris'] = surf1['nTris'] + surf2['nTris']
    surf['nVertices'] = surf1['nVertices'] + surf2['nVertices']
    surf['vertices'] = np.array(np.concatenate((surf1['vertices'],surf2['vertices'])))
    surf['tris'] = np.concatenate((surf1['tris'],surf2['tris']+surf1['nVertices']))
    return surf

@jit(nopython=True)
def alpha_from_beta(beta,aw,ap):
    """ Calculates the relative equi-distance alpha from relative equi-volume depth
    beta, using the local area of the bottom (white) and top (pial) vertex.
    """
    alpha = (2 * ap * beta - ap * beta**2 + aw * beta**2)/(ap + aw)
    return alpha

@jit(nopython=True)
def intermediate_plane(a1,b1,c1,a2,b2,c2,alpha_a,alpha_b,alpha_c):
    """ Calculates vertices of intermediate plane defined by relative depths for all three vertices independently
    """
    ma = a1 + alpha_a * (a2 - a1)
    mb = b1 + alpha_b * (b2 - b1)
    mc = c1 + alpha_c * (c2 - c1)
    return ma, mb, mc

@jit(nopython=True)
def equidist_plane(a1,b1,c1,a2,b2,c2,alpha):
    return intermediate_plane(a1,b1,c1,a2,b2,c2,alpha,alpha,alpha)

@jit(nopython=True)
def equivol_plane(a1,b1,c1,a2,b2,c2,beta,aw_a,aw_b,aw_c,ap_a,ap_b,ap_c):
    alpha_a = alpha_from_beta(beta,aw_a,ap_a);
    alpha_b = alpha_from_beta(beta,aw_b,ap_b);
    alpha_c = alpha_from_beta(beta,aw_c,ap_c);
    return intermediate_plane(a1,b1,c1,a2,b2,c2,alpha_a,alpha_b,alpha_c)

@jit(nopython=True)
def same_side(p1,p2,a,b):
    cp1 = np.cross(b-a, p1-a)
    cp2 = np.cross(b-a, p2-a)
    if np.dot(cp1, cp2) >= 0:
        return True
    else:
        return False

@jit(nopython=True)
def test_point_inside_triangle(p,a,b,c):
    if same_side(p,a,b,c) and same_side(p,b,a,c) and same_side(p,c,a,b):
        return True
    else:
        return False

@jit(nopython=True)
def plane_side(p,ma,mb,mc):
    return np.dot(p-ma,np.cross(mb-ma,mc-ma))

def find_depth_through_point(p,depth_plane_function):
    def f(x):
        ma, mb, mc = depth_plane_function(x)
        return plane_side(p,ma,mb,mc)
    x_min = minimize_scalar(f,bounds=[0,1],method='bounded').x
    x_max = minimize_scalar(lambda x: -f(x),bounds=[0,1],method='bounded').x
    x_grid = [0, min(x_min,x_max), max(x_min,x_max), 1]
    f_x_grid = [f(x_grid[0]), f(x_grid[1]), f(x_grid[2]), f(x_grid[3])]

    for i in [1,0,2]:
        if np.sign(f_x_grid[i])!=np.sign(f_x_grid[i+1]):
            x_zero = root_scalar(f,
                               method='brentq',
                               bracket=[x_grid[i],x_grid[i+1]]).root
            ma, mb, mc = depth_plane_function(x_zero)
            if test_point_inside_triangle(p,ma,mb,mc):
                return x_zero
    return None
    
@jit(nopython=True)
def numba_calc_bounding_prisms(faces,vertices_white,vertices_pial,n_x,n_y,n_z):
    bounding_prisms = np.zeros((n_x,n_y,n_z,1),dtype=np.int64)
    for i_prism, prism_vertices_indcs in enumerate(faces):
        prism_points = np.concatenate(
            (vertices_white[prism_vertices_indcs],
             vertices_pial[prism_vertices_indcs]))
        max_x_idx = int(np.round(np.max(prism_points[:,0])))
        min_x_idx = int(np.round(np.min(prism_points[:,0])))
        max_y_idx = int(np.round(np.max(prism_points[:,1])))
        min_y_idx = int(np.round(np.min(prism_points[:,1])))
        max_z_idx = int(np.round(np.max(prism_points[:,2])))
        min_z_idx = int(np.round(np.min(prism_points[:,2])))

        x_idx_range = np.arange(min_x_idx,max_x_idx+1)
        y_idx_range = np.arange(min_y_idx,max_y_idx+1)
        z_idx_range = np.arange(min_z_idx,max_z_idx+1)

        # remove those outside of the slab
        x_idx_range = x_idx_range[(x_idx_range>=0) & (x_idx_range<n_x)]
        y_idx_range = y_idx_range[(y_idx_range>=0) & (y_idx_range<n_y)]
        z_idx_range = z_idx_range[(z_idx_range>=0) & (z_idx_range<n_z)]

        for x_idx in x_idx_range:
            for y_idx in y_idx_range:
                for z_idx in z_idx_range:
                    available_indcs = np.argwhere(bounding_prisms[x_idx,y_idx,z_idx,:]==0)
                    if available_indcs.size==0:
                        old_length=int(bounding_prisms.shape[3])
                        bounding_prisms = np.concatenate((bounding_prisms,
                                                          np.zeros((n_x,n_y,n_z,old_length),
                                                                   dtype=np.int64)),
                                                         axis=3)
                        next_available_idx=old_length
                    else:
                        next_available_idx = int(available_indcs[0][0])                      
                    bounding_prisms[x_idx,y_idx,z_idx,next_available_idx] = int(i_prism)
    return bounding_prisms


def calc_depth_from_surfaces_voxelidx(voxelidx,faces,vertices_white,vertices_pial,
                                      area_white, area_pial,bounding_prisms,method):
    x_idx, y_idx, z_idx = voxelidx
    p = np.array(voxelidx)
    prism_idx=0
    while ((bounding_prisms[x_idx,y_idx,z_idx,prism_idx] != 0)
           and prism_idx<bounding_prisms.shape[3]):
        canditate_prism_idx = bounding_prisms[x_idx,y_idx,z_idx,prism_idx]
        a1 = vertices_white[faces[canditate_prism_idx,0],:]
        b1 = vertices_white[faces[canditate_prism_idx,1],:]
        c1 = vertices_white[faces[canditate_prism_idx,2],:]
        a2 = vertices_pial[faces[canditate_prism_idx,0],:]
        b2 = vertices_pial[faces[canditate_prism_idx,1],:]
        c2 = vertices_pial[faces[canditate_prism_idx,2],:]

        if method=='equivol':
            aw_a = area_white[faces[canditate_prism_idx,0]]
            aw_b = area_white[faces[canditate_prism_idx,1]]
            aw_c = area_white[faces[canditate_prism_idx,2]]
            ap_a = area_pial[faces[canditate_prism_idx,0]]
            ap_b = area_pial[faces[canditate_prism_idx,1]]
            ap_c = area_pial[faces[canditate_prism_idx,2]]            

            depth_plane_function = \
                lambda beta: equivol_plane(a1,b1,c1,a2,b2,c2,
                                           beta,
                                           aw_a,aw_b,aw_c,ap_a,ap_b,ap_c)
        elif method=='equidist':
            depth_plane_function = \
                lambda alpha: equidist_plane(a1,b1,c1,a2,b2,c2,alpha)
            
        d = find_depth_through_point(p,depth_plane_function)

        if d is not None:
            return d, float(canditate_prism_idx)

        prism_idx = prism_idx + 1
    return np.nan, np.nan

@jit(nopython=True)
def numba_assign_calc_depth_results(results,roi,grid_size):
    depths = np.empty(grid_size)
    columns = np.empty(grid_size)
    depths[:] = np.nan
    columns[:] = np.nan
    for roi_idx,result in zip(roi,results):
        depths[roi_idx[0],roi_idx[1],roi_idx[2]]= result[0]
        columns[roi_idx[0],roi_idx[1],roi_idx[2]]= result[1]
    return depths, columns
                                   
def calc_depth_from_surfaces_on_grid(surf_white, area_white, surf_pial, area_pial,
                                     n_x, n_y, n_z, method, n_jobs=32):
    faces = surf_pial['tris']
    vertices_white = surf_white['vertices']
    vertices_pial = surf_pial['vertices']

    bounding_prisms = numba_calc_bounding_prisms(faces,
                                                 vertices_white,vertices_pial,
                                                 n_x, n_y, n_z)

    roi = np.argwhere(np.any(bounding_prisms,axis=3))
    
    with parallel_backend('loky', inner_max_num_threads=1):
        results = List(
            Parallel(n_jobs=n_jobs)(
                delayed(
                    lambda roi_idx:
                    calc_depth_from_surfaces_voxelidx(
                        roi_idx,faces,vertices_white,vertices_pial,
                        area_white, area_pial,bounding_prisms,method))(roi_idx)
                for roi_idx in roi))

    depths, columns = numba_assign_calc_depth_results(
        results,roi,bounding_prisms.shape[:3])        

    return depths, columns


def process_dc_voxeldepth_from_surfaces(surf_white_lh_file,area_white_lh_file,
                                        surf_pial_lh_file,area_pial_lh_file,
                                        surf_white_rh_file,area_white_rh_file,
                                        surf_pial_rh_file,area_pial_rh_file,
                                        volume_file,
                                        depths_fname,columns_fname,
                                        method='equivol',
                                        upsample_factor=None,n_jobs=32):

    volume = nib.load(volume_file)
    n_x, n_y, n_z = volume.shape[:3]
    voxel_to_scanner = volume.affine

    # optionally use upsampled grid
    if upsample_factor is not None:
        n_x = int(np.floor(n_x * upsample_factor))
        n_y = int(np.floor(n_y * upsample_factor))
        n_z = int(np.floor(n_z * upsample_factor))
        upsampled_to_voxel = np.matrix(
            [[1/upsample_factor,0                ,0                ,-0.5 + 1/(2*upsample_factor)],
             [0                ,1/upsample_factor,0                ,-0.5 + 1/(2*upsample_factor)],
             [0                ,0                ,1/upsample_factor,-0.5 + 1/(2*upsample_factor)],
             [0                ,0                ,0                ,1]]);
        grid_to_scanner = voxel_to_scanner @ upsampled_to_voxel
    else:
        grid_to_scanner = voxel_to_scanner
        
    # load all surfaces in voxel space
    surf_white_lh = load_fs_surf_in_grid(surf_white_lh_file,grid_to_scanner)
    surf_pial_lh = load_fs_surf_in_grid(surf_pial_lh_file,grid_to_scanner)
    surf_white_rh = load_fs_surf_in_grid(surf_white_rh_file,grid_to_scanner)
    surf_pial_rh = load_fs_surf_in_grid(surf_pial_rh_file,grid_to_scanner)

    # load area files
    area_white_lh = nib.freesurfer.read_morph_data(area_white_lh_file)
    area_pial_lh = nib.freesurfer.read_morph_data(area_pial_lh_file)
    area_white_rh = nib.freesurfer.read_morph_data(area_white_rh_file)
    area_pial_rh = nib.freesurfer.read_morph_data(area_pial_rh_file)    

    # merge hemispheres
    surf_white = merge_surfaces(surf_white_lh,surf_white_rh)
    surf_pial = merge_surfaces(surf_pial_lh,surf_pial_rh)
    area_white = np.concatenate((area_white_lh,area_white_rh),axis=0)
    area_pial = np.concatenate((area_pial_lh,area_pial_rh),axis=0)        

    # calc voxel depths
    depths, columns = calc_depth_from_surfaces_on_grid(surf_white, area_white,
                                                       surf_pial, area_pial,
                                                       n_x, n_y, n_z,
                                                       method,n_jobs)

    xform = grid_to_scanner
    nii_depths = nib.nifti1.Nifti1Image(depths, xform)
    nii_columns = nib.nifti1.Nifti1Image(columns, xform)

    nib.save(nii_depths,depths_fname)
    nib.save(nii_columns,columns_fname)
    
    return nii_depths, nii_columns
