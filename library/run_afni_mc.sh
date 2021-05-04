#!/bin/bash

fNameBase=$(remove_ext $1)

nVols=$(fslinfo ${fNameBase} | grep ^dim4 | awk '{print $2}')

if [ -z "$2" ]
then
    base=$(expr ${nVols} / 2)
else
    base=$2
fi

3dvolreg -prefix ${fNameBase}_mc.nii \
         -Fourier \
         -float \
         -base $base \
         -dfile ${fNameBase}_mc.par \
         -maxdisp1D ${fNameBase}_mc_maxdisp \
         ${fNameBase}.nii

awk -F' ' '{s=$2;$1=$3;$2=-$4;$3=s;t=$5;$4=$6;$5=-$7;$6=t;$7=$8;$8=$9;$9=""}1' OFS=' ' \
    ${fNameBase}_mc.par >  ${fNameBase}_mc_reordered.par
fsl_tsplot -i ${fNameBase}_mc_reordered.par -t 'AFNI estimated rotations (radians)' \
           -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o ${fNameBase}_rot.png 
fsl_tsplot -i ${fNameBase}_mc_reordered.par -t 'AFNI estimated translations (mm)' \
           -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o ${fNameBase}_trans.png
fsl_tsplot -i ${fNameBase}_mc_maxdisp,${fNameBase}_mc_maxdisp_delt -t 'AFNI estimated max displacement (mm)' \
           -u 1 -w 640 -h 144 -a absolute,relative -o ${fNameBase}_disp.png
