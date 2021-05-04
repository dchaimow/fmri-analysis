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
boldcorrect.sh func_go-nogo ${shiftFraction}

cd ${curDir}
