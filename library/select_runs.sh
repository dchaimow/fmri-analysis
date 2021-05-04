#!/bin/bash
#
# select_runs.sh <runlist_file> <idx1> <idx2> ...
#
# - returns a selected list of files from runlist_file
# - lines idx1, idx2, ... are selected

select_runs_add-ext.sh $1 "" ${@:2}
