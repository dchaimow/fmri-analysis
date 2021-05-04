#!/bin/bash

subj=$1

dataDir=~/Data/FinnEtAl2019/ds002076

TR=4.012
shiftFraction=0.5



# import
importruns_vaso-split_reverse.sh \
    func ${dataDir}/${subj}/func/${subj}_task-layer_run-?_bold.nii.gz

# find out which tasks
find_task-runs.sh \
    func  ${dataDir}/${subj}/func/${subj}_task-layer_run-?_bold.nii.gz

# motioncorrect
motioncorrect_vaso.sh $(< func_runs_basenames.txt)

exit

# TODO: make sure to use the motion corrected runs
# avgruns according to conditions
avgruns func_alpha-rem_nulled.nii \
        $(select_runs_add-ext.sh func_runs_basename.txt nulled.nii \
                                 $(< func_alpha-rem_task-runs.txt))
avgruns func_alpha-rem_notnulled.nii \
        $(select_runs_add-ext.sh func_runs_basename.txt notnulled.nii \
                                 $(< func_alpha-rem_task-runs.txt))
avgruns func_go-nogo_nulled.nii \
        $(select_runs_add-ext.sh func_runs_basename.txt nulled.nii \
                                 $(< func_go-nogo_task-runs.txt))
avgruns func_go-nogo_notnulled.nii \
        $(select_runs_add-ext.sh func_runs_basename.txt notnulled.nii \
                                 $(< func_go-nogo_task-runs.txt))

# average all
avgruns func_all_nulled.nii func?_mc
        

# calct1

# bold correction

