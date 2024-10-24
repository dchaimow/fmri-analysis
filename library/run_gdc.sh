#!/bin/bash

file_in=$1
coeff_path=${2:-/data/hu_dchaimow/owncloud/pfc-layers/sequences/coeff_SC72CD.grad}

export FSLOUTPUTTYPE=NIFTI

curDir=$(pwd)
tmpdir=$(mktemp -d)
imcp ${file_in} ${tmpdir}/infile
cp ${coeff_path} ${tmpdir}/coeff.grad
cd ${tmpdir}

gradient_unwarp.py $(imglob -extension infile) corrected_tri.nii siemens -g coeff.grad -n
    
convertwarp --abs --ref=infile \
            --warp1=fullWarp_abs.nii.gz \
            --relout \
            --out=warpfield \
            --jacobian=jac

fslmaths jac -Tmean jac

applywarp --rel --interp=spline \
          -i infile \
          -r infile \
          -w warpfield \
          -o corrected \
          --paddingsize=1

# remove negative values from spline interpolation
fslmaths corrected -thr 0 corrected

cd ${curDir}

imcp ${tmpdir}/corrected $(remove_ext ${file_in})_gdc
imcp ${tmpdir}/warpfield $(remove_ext ${file_in})_gdc_warpfield
imcp ${tmpdir}/jac $(remove_ext ${file_in})_gdc_jac

rm -rf ${tmpDir}
