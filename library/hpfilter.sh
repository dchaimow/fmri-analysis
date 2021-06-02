#!/bin/bash
#
# hpfilter.sh <fName> <TR> <hpCutOff>
#
# - applies a temporal high-pass filter to fName
# - hpCutOff is the cut-off cycle (longest) in seconds 

fName=$1
TR=$2
hpCutOff=$3

hpSigma=$(bc -l <<< "(${hpCutOff}/2)/${TR}")
fslmaths $fName -Tmean tempMean
fslmaths $fName -bptf ${hpSigma} -1 -add tempMean $(remove_ext ${fName})_hpfilt
imrm tempMean


