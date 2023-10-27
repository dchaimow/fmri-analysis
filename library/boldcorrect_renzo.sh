#!/bin/bash
# Renzo's bold correction method (taken from afni_VASO_eval_SPM.sh)
#(adapted from a script written by Renzo Huber)
fBaseName=$1

# The first vaso volume is first nulled volume divided by the 2nd BOLD volume
3dcalc -prefix tmp_VASO_vol1.nii \
       -a      ${fBaseName}_notnulled.nii'[1]' \
       -b      ${fBaseName}_nulled.nii'[0]' \
       -expr 'b/a' -overwrite

# Calculate all VASO volumes after the first one
# -a goes from the 2nd BOLD volume to the 2nd-to-last BOLD volume
# -b goes from the 3rd BOLD volume to the last BOLD volume
# -c goes from the 2nd Nulled volume to the last Nulled volume
NumVol=`3dinfo -nv ${fBaseName}_notnulled.nii`
3dcalc -prefix tmp_VASO_othervols.nii \
       -a      ${fBaseName}_notnulled.nii'[0..'`expr $NumVol - 2`']' \
       -b      ${fBaseName}_notnulled.nii'[1..$]' \
       -c      ${fBaseName}_nulled.nii'[1..$]' \
       -expr 'c*2/(a+b)' -overwrite

# concatinate the first VASO volume with the rest of the sequence
3dTcat -overwrite -prefix ${fBaseName}_rvaso.nii \
       tmp_VASO_vol1.nii tmp_VASO_othervols.nii

# Remove the temporary separate files for the first VASO volume and the rest of the VASO volumes
rm tmp_VASO_vol1.nii tmp_VASO_othervols.nii

