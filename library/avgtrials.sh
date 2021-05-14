#!/bin/bash
#
# avgtrials.sh <filename> <avgname> <trial_duration> <TR> <trial_onsets ...>
#
# - event related trial averaging
# - uses 3ddeconvolve to estimate a finite impules model (TENT)

filename=$1
avgname=$2
trial_duration=$3
TR=$4
trial_onsets="${@:5}"

N=$(echo "scale=0;${trial_duration}/${TR}" | bc -l)
a=0
b=$(echo "${TR} * (${N}-1)" |  bc -l)
echo $b

3ddeconvolve -input ${filename} \
             -num_stimts 1 \
             -stim_times 1 "1D: ${trial_onsets}" "TENT($a,$b,$N)" \
             -polort a \
             -TR_times $TR \
             -x1D Model \
             -overwrite \
             -bout

# add mean estimate to impulse responses
NumVol=$(3dinfo -nv Decon+orig)
3dcalc -a Decon+orig'[1]' -b Decon+orig'['`expr $NumVol - $N`'..$]' \
       -expr 'a+b' -prefix ${avgname} -overwrite

# convert from multiple bricks to temporal nifti
3dTcat ${avgname}'[0..$]' -prefix ${avgname} -overwrite 
