#!/bin/bash

subj=$1

dataDir=$2

TR=4.012
shiftFraction=0.5

curDir=$(pwd)

mkdir -p ${subj}
cd ${subj}

# import
importruns_vaso-split_reverse.sh \
    func ${TR} ${dataDir}/${subj}/func/${subj}_task-layer_run-?_bold.nii.gz

# find out which tasks
find_task-runs.sh \
    func ${dataDir}/${subj}/func/${subj}_task-layer_run-?_bold.nii.gz

# motioncorrect
motioncorrect_vaso.sh $(< func_runs_basenames.txt)

# avgruns according to conditions
avgruns.sh func_alpha-rem_nulled.nii \
        $(select_runs_add-ext.sh func_runs_basenames.txt _nulled_mc.nii \
                                 $(< func_alpha-rem_task-runs.txt))
avgruns.sh func_alpha-rem_notnulled.nii \
        $(select_runs_add-ext.sh func_runs_basenames.txt _notnulled_mc.nii \
                                 $(< func_alpha-rem_task-runs.txt))
avgruns.sh func_go-nogo_nulled.nii \
        $(select_runs_add-ext.sh func_runs_basenames.txt _nulled_mc.nii \
                                 $(< func_go-nogo_task-runs.txt))
avgruns.sh func_go-nogo_notnulled.nii \
        $(select_runs_add-ext.sh func_runs_basenames.txt _notnulled_mc.nii \
                                 $(< func_go-nogo_task-runs.txt))

# average all (only for t1 calculation)
avgruns.sh func_all_nulled.nii func?_nulled_mc.nii
avgruns.sh func_all_notnulled.nii func?_notnulled_mc.nii

# calct1
calct1.sh func_all

# bold correction
boldcorrect.sh func_alpha-rem ${shiftFraction}
mv func_alpha-rem_notnulled.nii func_alpha-rem_bold.nii
boldcorrect.sh func_go-nogo ${shiftFraction}
mv func_go-nogo_notnulled.nii func_go-nogo_bold.nii

# # trial averaging
# trialOrder="1 0 0 0 1 1 0 1 0 0 1 1 1 0 1 0 0 0 1 1"
# onsetDelay_s=4 # stim starts 4s after this 4s delay
# trialDuration_s=32
# onsetTimes_A=$(python ../../code/calc_onset-times.py \
#                       0 ${onsetDelay_s} ${trialDuration} ${trialOrder})
# onsetTimes_B=$(python ../../code/calc_onset-times.py \
#                       1 ${onsetDelay_s} ${trialDuration} ${trialOrder})
# trialAvgDuration=44
# trial_dt=4
# for modality in bold vaso
# do
#     avgtrials.sh func_go-nogo_${modality}.nii trials_go_${modality}.nii \
#              $trialAvgDuration $trial_dt \
#              $onsetTimes_A
#     avgtrials.sh func_go-nogo_${modality}.nii trials_nogo_${modality}.nii \
#              $trialAvgDuration $trial_dt \
#              $onsetTimes_B

#     avgtrials.sh func_alpha-rem_${modality}.nii trials_alpha_${modality}.nii \
#              $trialAvgDuration $trial_dt \
#              $onsetTimes_A
#     avgtrials.sh func_alpha-rem_${modality}.nii trials_rem_${modality}.nii \
#              $trialAvgDuration $trial_dt \
#              $onsetTimes_B
# done


# # normalization
# for condition in go nogo alpha rem
# do
#     for modality in bold vaso
#     do
#         3dcalc -a trials_${condition}_${modality}.nii \
#                -b trials_${condition}_${modality}.nii'[0]' \
#                -expr '100*a/b - 100' \
#                -prefix trials_rchg_${condition}_${modality}.nii -overwrite
#     done
# done

# label=0002
# for layer in superficial deeper
# do
#     3dcalc -a lh.zstat1_smooth_clusters.label-${label}_in-func.nii \
#            -b ${layer}.nii \
#            -expr 'and(a,b)' \
#            -prefix layermask.nii -overwrite
#     for condition in go nogo alpha rem
#     do
#         for modality in bold vaso
#         do
#             3dmaskave -quiet \
#                       -mask layermask.nii \
#                       trials_rchg_${condition}_${modality}.nii \
#                       > ${modality}_${layer}_${condition}.1D
#         done
#     done
# done

# # visualize
# python3 ../../code/visualize_single-subject.py 

cd ${curDir}
