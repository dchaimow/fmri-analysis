#!/bin/bash
# bold correction taking asymmetric readout timing into account by doing linear interpolations
# based on Renzo's bold correction method (taken from afni_VASO_eval_SPM.sh)
# extends it by using readout timing dependent weights
# (adapted from a script written by Renzo Huber)

# boldcorrect_lin.sh <basename> <shiftFraction>
#
# - does bold correction of VASO by dividing nulled volumes by time shifted notnulled volumes
# - shiftFraction is shift relative to TR (usually a positive number)

fBaseName=$1
shiftFraction=$2

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
w1=$(bc -l <<< "1-${shiftFraction}")
w2=${shiftFraction}
3dcalc -prefix tmp_VASO_othervols.nii \
       -a      ${fBaseName}_notnulled.nii'[0..'`expr $NumVol - 2`']' \
       -b      ${fBaseName}_notnulled.nii'[1..$]' \
       -c      ${fBaseName}_nulled.nii'[1..$]' \
       -expr "c/(${w1}*a+${w2}*b)" -overwrite

# concatinate the first VASO volume with the rest of the sequence
3dTcat -overwrite -prefix ${fBaseName}_vaso.nii \
       tmp_VASO_vol1.nii tmp_VASO_othervols.nii

# Remove the temporary separate files for the first VASO volume and the rest of the VASO volumes
rm tmp_VASO_vol1.nii tmp_VASO_othervols.nii

