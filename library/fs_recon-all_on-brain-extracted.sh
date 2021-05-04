#!/bin/bash
#
# fs_recon-all_on-brain-extracted.sh <anatFileName> <subjectsDir> <subject>
#
# - runs freesurfer recon-all on high-res data that has already been brain extracted

anatFileName=$1
subjectsDir=$2
subject=$3

echo "mris_inflate -n 100" > expert.opts

recon-all -i ${anatFileName} \
          -hires \
          -autorecon1 \
          -noskullstrip \
          -sd ${subjectsDir} \
          -s ${subject} \
          -parallel    

ln -s \
   ${subjectsDir}/${subject}/mri/T1.mgz \
   ${subjectsDir}/${subject}/mri/brainmask.mgz

ln -s \
   ${subjectsDir}/${subject}/mri/T1.mgz \
   ${subjectsDir}/${subject}/mri/brainmask.auto.mgz  

recon-all -hires \
          -autorecon2 \
          -autorecon3 \
          -sd ${subjectsDir} \
          -s ${subject} \
          -expert expert.opts \
          -xopts-overwrite \
          -parallel   
