#!/bin/bash
#
# prepare_fieldmap <fmap_e1> <fmap_e2> <fmap_e2_ph> <anat_mp2rage_inv2>
#
# Prepares fieldmap.nii for distortion correction using fsl fugue.

fieldmap_e1=$1
fieldmap_e2=$2
fieldmap_ph=$3
mp2rage_inv2=$4

# register mp2rage INV2 to fieldmap magnitude of echo 1 (save transformation matrix only)
flirt -in ${mp2rage_inv2} \
      -ref ${fieldmap_e1} \
      -omat anat_to_fmap-e1.mat

# apply brain extraction to PD
bet ${mp2rage_inv2} anat_mp2rage-inv2_bet

# apply registration transformation to brain extracted PD
flirt -in anat_mp2rage-inv2_bet \
      -ref ${fieldmap_e1} \
      -applyxfm -init anat_to_fmap-e1.mat -out anat_mp2rage-inv2_bet_in-fmap   

# erode resulting brainmask
fslmaths anat_mp2rage-inv2_bet_in-fmap \
         -ero anat_mp2rage-inv2_bet_in-fmap_ero

# apply eroded brain mask to fieldmap magnitude of echo 1
fslmaths ${fieldmap_e1} \
         -mas anat_mp2rage-inv2_bet_in-fmap_ero \
         $(basename $(remove_ext ${fieldmap_e1})_bet)

# set deltaTE
TE_1=$(jq < $(remove_ext ${fieldmap_e1}).json '.EchoTime')
TE_2=$(jq < $(remove_ext ${fieldmap_e2}).json '.EchoTime')
deltaTE=$(echo "1000*($TE_2 - $TE_1)" | bc)

# prepare fieldmap
fsl_prepare_fieldmap SIEMENS \
    ${fieldmap_ph} \
    $(basename $(remove_ext ${fieldmap_e1}))_bet \
    fieldmap.nii \
    $deltaTE
