#!/usr/bin/env python

from skimage import morphology as sim
from scipy.ndimage.morphology import binary_fill_holes
from nipype.interfaces.freesurfer.preprocess import MRIConvert
from nipype.interfaces.fsl.maths import ApplyMask

import nibabel as nb
import numpy as np
import sys
import os

#
# anat_brain-extract_using-fs-refine.py <anat_fname> <fs_dir> <anatbrain_fname> <mask_fname>#

# refine_aseg() and grow_mask() are from fmriprep

def refine_aseg(aseg, ball_size=4):
    """
    Refine the ``aseg.mgz`` mask of Freesurfer.
    First step to reconcile ANTs' and FreeSurfer's brain masks.
    Here, the ``aseg.mgz`` mask from FreeSurfer is refined in two
    steps, using binary morphological operations:

    1. With a binary closing operation the sulci are included
       into the mask. This results in a smoother brain mask
       that does not exclude deep, wide sulci.

    2. Fill any holes (typically, there could be a hole next to
       the pineal gland and the corpora quadrigemina if the great
       cerebral brain is segmented out).
    """
    from skimage import morphology as sim
    from scipy.ndimage.morphology import binary_fill_holes

    # Read aseg data
    bmask = aseg.copy()
    bmask[bmask > 0] = 1
    bmask = bmask.astype(np.uint8)

    # Morphological operations
    selem = sim.ball(ball_size)
    newmask = sim.binary_closing(bmask, selem)
    newmask = binary_fill_holes(newmask.astype(np.uint8), selem).astype(np.uint8)

    return newmask.astype(np.uint8)

def grow_mask(anat, aseg, ants_segs=None, ww=7, zval=2.0, bw=4):
    """
    Grow mask including pixels that have a high likelihood.

    GM tissue parameters are sampled in image patches of ``ww`` size.
    This is inspired on mindboggle's solution to the problem:
    https://github.com/nipy/mindboggle/blob/master/mindboggle/guts/segment.py#L1660

    """
    from skimage import morphology as sim

    selem = sim.ball(bw)

    if ants_segs is None:
        ants_segs = np.zeros_like(aseg, dtype=np.uint8)

    aseg[aseg == 42] = 3  # Collapse both hemispheres
    gm = anat.copy()
    gm[aseg != 3] = 0

    refined = refine_aseg(aseg)
    newrefmask = sim.binary_dilation(refined, selem) - refined
    indices = np.argwhere(newrefmask > 0)
    for pixel in indices:
        # When ATROPOS identified the pixel as GM, set and carry on
        if ants_segs[tuple(pixel)] == 2:
            refined[tuple(pixel)] = 1
            continue

        window = gm[
            pixel[0] - ww:pixel[0] + ww,
            pixel[1] - ww:pixel[1] + ww,
            pixel[2] - ww:pixel[2] + ww,
        ]
        if np.any(window > 0):
            mu = window[window > 0].mean()
            sigma = max(window[window > 0].std(), 1.0e-5)
            zstat = abs(anat[tuple(pixel)] - mu) / sigma
            refined[tuple(pixel)] = int(zstat < zval)

    refined = sim.binary_opening(refined, selem)
    return refined

if __name__ == "__main__":
    if len(sys.argv)!=5:
        print("Usage: anat_brain-extract_using-fs-refine.py <anat_fname> <fs_dir> <anatbrain_out_fname> <mask_out_fname>");
        sys.exit()
    
    anat_fname = sys.argv[1]
    fs_dir = sys.argv[2]
    anatbrain_fname = sys.argv[3]
    mask_fname = sys.argv[4]

    # convert aseg image from freesurfer stream
    result = MRIConvert(in_file=os.path.join(fs_dir,'mri','aseg.mgz'),
                        out_orientation='RAS',
                        resample_type='nearest',
                        reslice_like=anat_fname).run()

    # load aseg and t1 anat
    anatnii = nb.load(anat_fname)
    asegnii= nb.load(result.outputs.out_file)
    
    # calculate refined brain mask
    masknii = nb.Nifti1Image(
        grow_mask(
            anatnii.get_fdata(dtype="float32"),
            np.asanyarray(asegnii.dataobj).astype("int16")),
        anatnii.affine,
        anatnii.header)
    masknii.set_data_dtype(np.uint8)

    masknii.to_filename(mask_fname)
    ApplyMask(in_file=anat_fname, mask_file=mask_fname,
              out_file=anatbrain_fname).run()
