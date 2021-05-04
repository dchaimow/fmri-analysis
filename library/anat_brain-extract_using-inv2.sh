#!/bin/bash
#
# anat_brain-extract_using-inv2.sh <anat> <anat_mp2rage_inv2>
#
# - brain extraction of anatomy (MP2RAGE UNI or processed) by running
#   bet on inv2 from the same acquisition and applying the resulting mask

anat=$1
anat_mp2rage_inv2=$2

bet $anat_mp2rage_inv2 $(basename $(remove_ext $anat_mp2rage_inv2))_brain

fslmaths $anat \
         -mas $(basename $(remove_ext $anat_mp2rage_inv2))_brain \
         $(basename $(remove_ext $anat))_brain
