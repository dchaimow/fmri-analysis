#!/bin/bash
#
# importruns.sh <basename> <run1_file> <run2_file>
#
# -imports functional runs into current directory and names them as
#  basename1 basename2 ...
# - uses fslmaths for copying, setting output data type to float
# -writes list of imported file names into runs.txt

fBase=$1
inFileNames="${@:2}"

i=0
> runs.txt
for inFileName in $inFileNames
do
    i=$(( i+1 ))
    
    echo ${i}: $inFileName 
    fBaseName=${fBase}$i

    fslmaths ${inFileName} ${fBaseName}.nii -odt float
    echo ${fBaseName}.nii >> runs.txt
done
