from nilearn import plotting as nl_plotting
from nilearn._utils import check_niimg
import numpy as np
import nibabel as nib
import copy


def img_to_voxel_space(img):
    if img is not None:
        img = check_niimg(img, dtype="auto")
        img_voxel_space = copy.deepcopy(img)
        xform = np.eye(4)
        if np.linalg.det(img.affine) < 0:
            xform[0, 0] = -1
        img_voxel_space.set_qform(xform)
        img_voxel_space.set_sform(xform)
        return img_voxel_space
    else:
        return None


def plot_anat(img, cut_coords=None, axes=None, dim=0):
    img = img_to_voxel_space(img)
    return nl_plotting.plot_anat(
        img, display_mode="z", cut_coords=cut_coords, axes=axes, dim=dim
    )


def plot_epi(img, cut_coords=None, axes=None):
    img = img_to_voxel_space(img)
    return nl_plotting.plot_epi(img, display_mode="z", cut_coords=cut_coords, axes=axes)


def plot_stat_map(stat_map_img, bg_img, threshold=1e-06, axes=None):
    stat_map_img = img_to_voxel_space(stat_map_img)
    bg_img = img_to_voxel_space(bg_img)
    return nl_plotting.plot_stat_map(
        stat_map_img=stat_map_img,
        bg_img=bg_img,
        display_mode="z",
        threshold=threshold,
        axes=axes,
    )


def plot_roi(
    roi_img,
    bg_img,
    display_mode="z",
    cut_coords=None,
    cmap=None,
    dim=0,
    axes=None,
    alpha=0.5,
    annotate=True,
    title=None,
    radiological=False,
):
    roi_img = img_to_voxel_space(roi_img)
    bg_img = img_to_voxel_space(bg_img)
    return nl_plotting.plot_roi(
        roi_img=roi_img,
        bg_img=bg_img,
        display_mode=display_mode,
        cmap=cmap,
        cut_coords=cut_coords,
        dim=dim,
        axes=axes,
        alpha=alpha,
        annotate=annotate,
        title=title,
    )


#   nl_plotting.plot_roi(roi_img=roi_img,bg_img=bg_img,display_mode='z',colorbar=True,cmap='Paired')


def plot_labels(
    roi_img, bg_img, display_mode="z", cut_coords=None, cmap=None, axes=None
):
    roi_img = img_to_voxel_space(roi_img)
    bg_img = img_to_voxel_space(bg_img)
    img = roi_img.get_fdata()
    labels = np.unique(img)
    img_new_labels = np.zeros(img.shape, dtype=np.int)
    for i, label in enumerate(labels):
        if i != 0:
            img_new_labels[img == label] = i
            print(i, label)
    roi = nib.Nifti2Image(img_new_labels, roi_img.affine)
    return nl_plotting.plot_roi(
        roi_img=roi,
        bg_img=bg_img,
        display_mode=display_mode,
        cmap=cmap,
        cut_coords=cut_coords,
        axes=axes,
    )


#   nl_plotting.plot_roi(roi_img=roi_img,bg_img=bg_img,display_mode='z',colorbar=True,cmap='Paired')
