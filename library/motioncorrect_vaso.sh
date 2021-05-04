#!/bin/bash
#
# motioncorrect_vaso.sh <run1_basename> <run2_basename> ...
#
# - runs motion correction on a list of vaso runs
# - registering them all to a common robust volume
# - uses afni and depends on run_afni_mc.sh
# - writes output as <run1>_[nulled|notnulled]_mc.nii <run2>_[nulled|notnulled]_mc.nii ...


fBaseNames="$@"

# find best volume to register to
nVols=()
> all_nulled_outcount.1D
> all_notnulled_outcount.1D
for fBaseName in $fBaseNames
do
    3dToutcount -automask -fraction -polort 5 -legendre \
                ${fBaseName}_nulled.nii'[3..$]' >> all_nulled_outcount.1D
    3dToutcount -automask -fraction -polort 5 -legendre \
                ${fBaseName}_notnulled.nii'[3..$]' >> all_notnulled_outcount.1D 
    nVols+=($(expr $(3dinfo -nti ${fBaseName}_nulled.nii) +  1 - 3))
done

# add nulled and notnulled outlier counts
1deval -a  all_nulled_outcount.1D -b   all_notnulled_outcount.1D  -expr 'a + b' \
       > all_outcount_sum.1D

## get run number and TR index for minimum outlier volume
minindex=`3dTstat -argmin -prefix -  all_outcount_sum.1D\'`
ovals=(`1d_tool.py -set_run_lengths ${nVols[@]} -index_to_run_tr $minindex`)

## save run and TR indices for extraction of vr_base_min_outlier
minoutrun=${ovals[0]}
minouttr=$(expr ${ovals[1]} + 2)
echo "min outlier: run $minoutrun, TR $minouttr" | tee min_outlier.txt

# run motion correction
fBaseNameList=($fBaseNames)
mcbaseNulled=${fBaseNameList[$(expr ${minoutrun} - 1)]}_nulled.nii\[${minouttr}\]
mcbaseNotNulled=${fBaseNameList[$(expr ${minoutrun} - 1)]}_notnulled.nii\[${minouttr}\]
for fBaseName in $fBaseNames
do
    run_afni_mc.sh ${fBaseName}_nulled ${mcbaseNulled}
    run_afni_mc.sh ${fBaseName}_notnulled ${mcbaseNotNulled}
done
