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
3dTstat -cvarinvNOD -prefix ${fBaseName}_T1_raw.nii \
        -overwrite ${fBaseName}_combined.nii 
rm ${fBaseName}_combined.nii

3dTstat -mean -prefix ${fBaseName}_mean_notnulled.nii ${fBaseName}_notnulled.nii -overwrite
3dTstat -mean -prefix ${fBaseName}_mean_nulled.nii ${fBaseName}_notnulled.nii -overwrite

3drefit -space ORIG -view orig ${fBaseName}_mean_notnulled.nii
3drefit -space ORIG -view orig ${fBaseName}_mean_nulled.nii

LN_MP2RAGE_DNOISE -INV1 ${fBaseName}_mean_nulled.nii -INV2 ${fBaseName}_mean_nulled.nii \
                  -UNI ${fBaseName}_T1_raw.nii -beta 5 -output ${fBaseName}_T1_denoised.nii

N4BiasFieldCorrection -i ${fBaseName}_T1_denoised.nii -o ${fBaseName}_T1.nii
rm ${fBaseName}_T1_raw.nii
rm ${fBaseName}_T1_denoised.nii
