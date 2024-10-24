#!/usr/local/fsl/bin/bash
# Spatial Smoothing using SUSAN (edge preserving)
fBaseName=$(basename $(basename $1 .gz) .nii)
FWHMsmooth=$2

export FSLOUTPUTTYPE=NIFTI

# create a mask based on all all time points being above 10% of the range between the 2% and 98% percentile
prctiles=($(fslstats ${fBaseName} -p 2 -p 98))
threshold=$(bc -l <<< "${prctiles[0]} + (${prctiles[1]}-${prctiles[0]})/10")
fslmaths ${fBaseName} -thr ${threshold} -Tmin -bin ${fBaseName}_mask -odt char

# calculate 50% percentile of voxels inside mask
medianInMask=$(fslstats ${fBaseName} -k ${fBaseName}_mask -p 50)

# dilate mask
fslmaths ${fBaseName}_mask -dilF ${fBaseName}_mask

# apply mask to data
fslmaths ${fBaseName} -mas ${fBaseName}_mask ${fBaseName}_masked

# calculate mean func
fslmaths ${fBaseName}_masked -Tmean ${fBaseName}_mean

# SUSAN - edge preserving spatial smoothing
# bright threshold = 0.75 * (50% precentile above - 2% percentile above)
brightThreshold=$(bc -l <<< "0.75 * (${medianInMask}-${prctiles[1]})") 
sigmaSmooth=$(bc -l <<< "${FWHMsmooth}/2.355")
susan ${fBaseName} ${brightThreshold} ${sigmaSmooth} 3 1 1 ${fBaseName}_mean ${brightThreshold} ${fBaseName}_smooth

# (re)apply mask to smoothed data
fslmaths ${fBaseName}_smooth -mas ${fBaseName}_mask ${fBaseName}_smooth
