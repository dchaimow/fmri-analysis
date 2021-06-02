#!/bin/bash
#
# importruns_vaso.sh <basename> <TR> <run1_file> <run2_file>
#
# - imports already splitted functional vaso runs into current assuming nulled to be acquired first
# - sets the TR
# - overwrites 1st two volumes of each
# - writes list of imported base file names into runs_basenames.txt

fBase=$1
TR=$2
inFileBaseNames="${@:3}"

i=0
> ${fBase}_runs_basenames.txt
for inFileBaseName in $inFileBaseNames
do
    i=$(( i+1 ))
    
    echo ${i}: $inFileName 
    fBaseName=${fBase}$i
    
    # import files
    3dTcat -prefix ${fBaseName}_nulled.nii ${inFileBaseName}_nulled.nii
    3dTcat -prefix ${fBaseName}_notnulled.nii ${inFileBaseName}_notnulled.nii    
    
    # set TR
    3drefit -TR $TR ${fBaseName}_nulled.nii
    3drefit -TR $TR ${fBaseName}_notnulled.nii
    
    # overwrite 1st 2 volumes (each)
    3dTcat -overwrite -prefix ${fBaseName}_nulled.nii \
           ${fBaseName}_nulled.nii'[2..3]' \
           ${fBaseName}_nulled.nii'[2..$]'
    3dTcat -overwrite -prefix ${fBaseName}_notnulled.nii \
           ${fBaseName}_notnulled.nii'[2..3]' \
           ${fBaseName}_notnulled.nii'[2..$]'

    echo ${fBaseName} >> ${fBase}_runs_basenames.txt
done
