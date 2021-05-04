#!/bin/bash
#
# importruns_vaso-split_reverse.sh <basename> <TR> <run1_file> <run2_file>
#
# - imports functional vaso runs into current directory and splits them into nulled and
#   non-nulled, ASSUMING NOTNULLED TO BE FIRST!
# - sets the TR
# - overwrites 1st two volumes of each
# - writes list of imported base file names into runs_basenames.txt

fBase=$1
TR=$2
inFileNames="${@:3}"

i=0
> runs_basenames.txt
for inFileName in $inFileNames
do
    i=$(( i+1 ))
    
    echo ${i}: $inFileName 
    fBaseName=${fBase}$i
    
    # split into nulled and bold (not nulled)
    3dTcat -prefix ${fBaseName}_nulled.nii ${inFileName}'[1..$(2)]'
    3dTcat -prefix ${fBaseName}_notnulled.nii ${inFileName}'[2..$(2)]'
    
    # set TR
    3drefit -TR $TR ${fBaseName}_nulled.nii
    3drefit -TR $TR ${fBaseName}_notnulled.nii
    
    # overwrite 1st 2 volumes (each)
    3dTcat -overwrite -prefix ${fBaseName}_nulled.nii \
           ${fBaseName}_nulled.nii'[2..3]' \
           ${fBaseName}_nulled.nii'[2..$]'
    3dTcat -overwrite -prefix ${fBaseName}_bold.nii \
           ${fBaseName}_notnulled.nii'[2..3]' \
           ${fBaseName}_notnulled.nii'[2..$]'

    echo ${fBaseName} >> runs_basenames.txt
done
