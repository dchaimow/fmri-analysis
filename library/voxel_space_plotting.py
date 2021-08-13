from nilearn import plotting as nl_plotting
from nilearn._utils import check_niimg
import numpy as np
import nibabel as nib
import copy

def img_to_voxel_space(img):
    if img is not None:
        img = check_niimg(img,dtype='auto')
        img_voxel_space = copy.deepcopy(img)        
        xform = np.eye(4)
        if np.linalg.det(img.affine)<0:
            xform[0,0] = -1
        img_voxel_space.set_qform(xform)
        img_voxel_space.set_sform(xform)
        return img_voxel_space
    else:
        return None
    
def plot_anat(img):
    img = img_to_voxel_space(img)
    nl_plotting.plot_anat(img,display_mode='z')

def plot_epi(img):
    img = img_to_voxel_space(img)
    nl_plotting.plot_epi(img,display_mode='z')

def plot_stat_map(stat_map_img,bg_img,threshold=1e-06):
    stat_map_img = img_to_voxel_space(stat_map_img)
    bg_img = img_to_voxel_space(bg_img)
    nl_plotting.plot_stat_map(stat_map_img=stat_map_img,bg_img=bg_img,display_mode='z',threshold=threshold)

def plot_roi(roi_img,bg_img):
    roi_img = img_to_voxel_space(roi_img)
    bg_img = img_to_voxel_space(bg_img)
    nl_plotting.plot_roi(roi_img=roi_img,bg_img=bg_img,display_mode='z')
#    nl_plotting.plot_roi(roi_img=roi_img,bg_img=bg_img,display_mode='z',colorbar=True,cmap='Paired')

