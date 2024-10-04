#!/usr/local/fsl/bin/bash

# glm analysis for a simple [rest activation] x repeat paradigm

basename=$1
T_rest=$2
T_act=$3

TR=$(3dinfo -tr ${basename}.nii)
NumVol=$(3dinfo -nv ${basename}.nii)
run_duration=$(bc -l <<< "${TR}*${NumVol}")
ITI=$(bc -l <<< "${T_rest}+${T_act}")

3dDeconvolve -input ${basename}_preproc_bold.nii \
             -polort A \
             -num_stimts 1 -stim_times 1 "1D: $(seq ${Ton} ${ITI} ${run_duration})" "BLOCK(${T_act},1)" \
             -CENSORTR 0-2 \
             -tout -bucket Decon.nii -overwrite
3dTcat -prefix ${basename}_tstat.nii Decon.nii'[2]' -overwrite

