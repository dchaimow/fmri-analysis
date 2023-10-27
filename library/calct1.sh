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

# clip 1 - make bias field correction less dependent on few extreme voxels
3dcalc -datum short -nscale -a ${fBaseName}_T1_denoised.nii \
       -expr "32767*max(0,min(1,a/(10*$(3dClipLevel ${fBaseName}_T1_denoised.nii))))" \
       -prefix ${fBaseName}_T1_denoised.nii -overwrite
                                  
# bias field correction
spm_bias-correct.py ${fBaseName}_T1_denoised.nii
mv m${fBaseName}_T1_denoised.nii ${fBaseName}_T1.nii
fslcpgeom ${fBaseName}_T1_denoised.nii ${fBaseName}_T1.nii # spm appears to change the affine slightly
#N4BiasFieldCorrection -i ${fBaseName}_T1_denoised.nii -o ${fBaseName}_T1.ni

# clip 2 - clip values that exceed expected distribution of intensities
3dcalc -datum byte -nscale -a ${fBaseName}_T1.nii \
       -expr "255*max(0,min(1,a/(4.2*$(3dClipLevel ${fBaseName}_T1.nii))))" \
       -prefix ${fBaseName}_T1.nii -overwrite

# calculate mask from notnulled
export FSLOUTPUTTYPE=NIFTI
3dAutomask -dilate 1 -prefix ${fBaseName}_mask.nii ${fBaseName}_notnulled.nii -overwrite
# apply mask to T1
fslmaths ${fBaseName}_T1.nii -mas  ${fBaseName}_mask.nii  ${fBaseName}_T1_brain.nii

rm ${fBaseName}_T1_raw.nii
rm ${fBaseName}_T1_denoised.nii
rm ${fBaseName}_T1_denoised_border_enhance.nii
rm ${fBaseName}_mask.nii
