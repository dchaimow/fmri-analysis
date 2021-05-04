#!/bin/bash
#
# select_runs_add-ext.sh <runlist_file> <ext> <idx1> <idx2> ...
#
# - returns a selected list of files from runlist_file with added extension ext
# - e.g. when file contains only basename and nulle/nonnulled needs to be added
# - lines idx1, idx2, ... are selected

runlistFile=$1
ext=$2
fileIdcs=${@:3}

fNames=($(< $runlistFile))

selected_fnames=''
for i in $fileIdcs
do
    ((l=i-1))
    selected_fnames+="${fNames[$l]}${ext} "
done
echo $selected_fnames
