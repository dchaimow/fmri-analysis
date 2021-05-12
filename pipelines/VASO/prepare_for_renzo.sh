#!/bin/bash

sessionDir=$1
inFileBase=$2

analysisDir=${3:-analysis_for-renzo}

TR=${4:-3.70202}
TR1=${5:-1.51440}
shiftFraction=$(bc -l <<< "${TR1}/${TR}")

export FSLOUTPUTTYPE=NIFTI

curDir=$(pwd)

mkdir -p ${sessionDir}/${analysisDir}
cd ${sessionDir}/${analysisDir}

# import raw
cp ${sessionDir}/func/${inFileBase}_nulled.nii func_nulled.nii
cp ${sessionDir}/func/${inFileBase}_notnulled.nii func_notnulled.nii
cp func_nulled.nii func-raw_nulled.nii
cp func_notnulled.nii func-raw_notnulled.nii

# set TR, overwrite 1st two volumes
3drefit -TR $TR func_nulled.nii
3drefit -TR $TR func_notnulled.nii

# overwrite 1st 2 volumes (each)
3dTcat -overwrite -prefix func_nulled.nii \
       func_nulled.nii'[2..3]' \
       func_nulled.nii'[2..$]'
3dTcat -overwrite -prefix func_notnulled.nii \
           func_notnulled.nii'[2..3]' \
           func_notnulled.nii'[2..$]'

# motioncorrect and overwrite
motioncorrect_vaso.sh func
mv func_nulled_mc.nii func_nulled.nii
mv func_notnulled_mc.nii func_notnulled.nii
rm -f *_mc_*
rm -f *.1D
rm -f *.png
rm -f min_outlier.txt

# bold correction
boldcorrect.sh func ${shiftFraction}
rm -f func_notnulled_tshift.nii

# calct1
calct1.sh func
rm -f func_mask.nii
rm -f func_T1_brain.nii

# calc QA metrics
LN_SKEW -input func_nulled.nii
LN_SKEW -input func_notnulled.nii
LN_SKEW -input func_vaso.nii

cd ${curDir}
