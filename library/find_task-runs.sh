#!/bin/bash


fBase=$1
inFileNames="${@:2}"

rm -f ${fBase}_*_task-runs.txt
i=0
for inFileName in $inFileNames
do
    i=$(( i+1 ))
    task=$(jq < $(remove_ext ${inFileName}).json '.TaskName')
    # remove quotes (first suffix then prefix):
    task="${task%\"}"
    task="${task#\"}"
    echo "$i " >> ${fBase}_${task}_task-runs.txt
done

