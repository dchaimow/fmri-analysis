#!/bin/bash
#
# register_fs-to-vasoT1_prepare.sh <vaso_T1_file> <fs_dir> [<itksnap_binary>]
#
# - converts FS T1 to nifti
# - starts ITK-SNAP in order to perform semi-automatic rigid-body registration in ITK-SNAP and save transformation matrix as initial_matrix.txt
# - after running this, register_fs-to-vasoT1.sh will run non-linear registration without further manual intervention

vaso_T1_file=$1
fs_dir=$2
itksnap_binary=${3:-itksnap}

mri_convert ${fs_dir}/mri/T1.mgz fs_T1.nii

echo "Please perform semi-automatic rigid-body registration in ITK-SNAP and save transformation matrix as initial_matrix.txt"
${itksnap_binary} -g ${vaso_T1_file} -o fs_T1.nii   
