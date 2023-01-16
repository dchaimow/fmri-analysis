#!/bin/bash
#
# register_fs-to-vasoT1_no-manual.sh <vaso_T1_file> <fs_dir>
#
# - converts FS T1 to nifti
# - if initial_matrix.txt does not exist: starts ITK-SNAP in order to perform semi-automatic rigid-body registration in ITK-SNAP and save transformation matrix as initial_matrix.txt
# - runs non-linear registration using ants

vaso_T1_file=$1
fs_dir=$2

mri_convert ${fs_dir}/mri/T1.mgz fs_T1.nii

antsRegistration \
    --verbose 1 \
    --dimensionality 3  \
    --float 0  \
    --collapse-output-transforms 1  \
    --interpolation BSpline[5] \
    --output [fs_to_func_,fs_to_func_Warped.nii,fs_to_func_InverseWarped.nii]  \
    --use-histogram-matching 0  \
    --winsorize-image-intensities [0.005,0.995]  \
    --transform Rigid[0.1]  \
    --metric MI[${vaso_T1_file},fs_T1.nii,1,32,Regular,0.25]  \
    --convergence [1000x500x250x100,1e-6,10]  \
    --shrink-factors 12x8x4x2  \
    --smoothing-sigmas 4x3x2x1vox  \
    -x mask.nii \
    --transform Affine[0.1]  \
    --metric MI[${vaso_T1_file},fs_T1.nii,1,32,Regular,0.25]  \
    --convergence [1000x500x250x100,1e-6,10]  \
    --shrink-factors 12x8x4x2  \
    --smoothing-sigmas 4x3x2x1vox  \
    -x mask.nii \
    --transform SyN[0.1,3,0]  \
    --metric CC[${vaso_T1_file},fs_T1.nii,1,4]  \
    --convergence [50x50x70x50x20,1e-6,10]  \
    --shrink-factors 10x6x4x2x1  \
    --smoothing-sigmas 5x3x2x1x0vox  \
    -x mask.nii

cp fs_to_func_Warped.nii fs_t1_in-func.nii
fslcpgeom ${vaso_T1_file} fs_t1_in-func.nii # correct for possible small affine changes
