#!/bin/bash
#
# register_fs-to-vasoT1_no-manual.sh <vaso_T1_file> <fs_dir>
#
# - converts FS T1 to nifti
# - uses init.txt
# - runs non-linear registration using ants

bold_file=$1
fs_dir=$2
cwd=$3

mri_convert ${fs_dir}/mri/brain.mgz fs_brain.nii

ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS

n4bold_file=$(remove_ext ${bold_file})_n4.nii
N4BiasFieldCorrection -i ${bold_file} -o ${n4bold_file}
bold_file=${n4bold_file}

bold_brain_file=$(remove_ext ${bold_file})_brain.nii
3dAutomask -apply_prefix ${bold_brain_file} ${bold_file}
bold_file=${bold_brain_file}

antsRegistration \
    --verbose 1 \
    --dimensionality 3 \
    --float 0 \
    --output [fs_to_func_,fs_to_func_Warped.nii,fs_to_func_InverseWarped.nii]  \
    --use-histogram-matching 0 \
    --winsorize-image-intensities [0.005,0.995] \
    --initial-moving-transform init.txt \
    --transform Rigid[0.1] \
    --metric MI[${bold_file},fs_brain.nii,1,32,Regular,0.25] \
    --convergence [1000x500x250,1e-6,10] \
    --shrink-factors 4x2x1 \
    --smoothing-sigmas 2x1x0vox \
    --transform Affine[0.1] \
    --metric MI[${bold_file},fs_brain.nii,1,32,Regular,0.25] \
    --convergence [1000x500x250x100,1e-6,10] \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    --transform SyN[0.1,3,0] \
    --metric CC[${bold_file},fs_brain.nii,1,4] \
    --convergence [100x100x100x100,1e-6,10] \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    -x mask.nii > antsRegistration.log 2>&1

cp fs_to_func_Warped.nii fs_t1_in-func.nii
fslcpgeom ${bold_file} fs_t1_in-func.nii # correct for possible small affine changes
