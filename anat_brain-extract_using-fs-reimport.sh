#!/bin/bash
#
# anat_brain-extract_using-fs-reimport.sh <anat> <fsdir>
#
# - brain extraction of processed anatomy reimporting brain extracted freesurfer
#   volume using coordinates of mp2rage and using it as a mask

anat=$1
fsdir=$2


mri_convert --out_orientation RAS -rt nearest --reslice_like ${anat} \
            ${fsdir}/mri/brainmask.mgz fs-brainmask_in-orig_brain.nii

fslmaths fs-brainmask_in-orig_brain.nii -abs -bin \
         brainmask.nii

fslmaths $anat -mas brainmask $(basename $(remove_ext ${anat}))_brain
