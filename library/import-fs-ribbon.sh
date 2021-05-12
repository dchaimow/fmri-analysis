#!/bin/bash
#
# import-fs-ribbon.sh <freesurferDir> <analysisDir> <vaso_t1>
#
# - imports ribbon file from freesurfer
# - tranaforms it to func space
# - makes values compatible with LAYNII

freesurferDir=$1
analysisDir=$2
vaso_t1=$3

mri_convert ${freesurferDir}/mri/ribbon.mgz ${analysisDir}/fs_ribbon.nii
antsApplyTransforms --interpolation BSpline[5] \
                    -d 3 \
                    -i ${analysisDir}/fs_ribbon.nii \
                    -r ${vaso_t1} \
                    -t ${analysisDir}/fs_to_func_1Warp.nii.gz \
                    -t ${analysisDir}/fs_to_func_0GenericAffine.mat \
                    -o ${analysisDir}/fs_ribbon_in-func.nii \
                    -n NearestNeighbor
fslmaths ${analysisDir}/fs_ribbon_in-func.nii -rem 39 -max 1 \
         ${analysisDir}/rim.nii
