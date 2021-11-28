#!/bin/bash
#
# calc1t.sh <basename>
#
# - calculate t1 contrast from VASO data (nulled and notnulled according to basename)
# - also calculates mask from notnulled data and generated an additional brain only T1 image

fBaseName=$1
NumVol=$(3dinfo -nv ${fBaseName}_nulled.nii)

# combine nulled and notnulled
3dTcat -prefix ${fBaseName}_combined.nii  \
       ${fBaseName}_nulled.nii'[3..'`expr $NumVol - 2`']' \
       ${fBaseName}_notnulled.nii'[3..'`expr $NumVol - 2`']'
3dTstat -cvarinvNOD -prefix ${fBaseName}_T1.nii \
        -overwrite ${fBaseName}_combined.nii 
rm ${fBaseName}_combined.nii

# calculate mask from notnulled
3dAutomask -dilate 1 -prefix ${fBaseName}_mask.nii ${fBaseName}_notnulled.nii
# applt mask to T1
fslmaths ${fBaseName}_T1.nii -mas  ${fBaseName}_mask.nii  ${fBaseName}_T1_brain.nii
