#!/bin/bash

# usage: find_max-roi_slice.sh <roi_file> <dimension>
# roi_file: the file containing the roi
# dimension: the dimension to slice along (x, y, or z)

roi_file=$1
dimension=$2

if [[ $dimension == "x" ]]; then
    dimn=dim1
elif [[ $dimension == "y" ]]; then
    dimn=dim2
elif [[ $dimension == "z" ]]; then
    dimn=dim3
else
    echo "Invalid dimension. Please use x, y, or z."
    exit 1
fi

cur_dir=$(pwd)
tmp_dir=$(mktemp -d)
cp $roi_file $tmp_dir
cd $tmp_dir

# get the number of slices in the roi file
number_of_slices=$(fslval ${roi_file} $dimn)

# loop through each slice, count the number of voxels and save in a list
max_slice=0
max_num_voxels=0
for i in $(seq 0 $((number_of_slices-1))); do
    if [[ $dimension == "x" ]]; then
        fslroi ${roi_file} roi_slice_${i}.nii.gz $i 1 0 -1 0 -1
    elif [[ $dimension == "y" ]]; then
        fslroi ${roi_file} roi_slice_${i}.nii.gz 0 -1 $i 1 0 -1
    elif [[ $dimension == "z" ]]; then
        fslroi ${roi_file} roi_slice_${i}.nii.gz 0 -1 0 -1 $i 1
    fi
    num_voxels=$(fslstats roi_slice_${i}.nii.gz -V | awk '{print $1}')
    if [[ $num_voxels -gt $max_num_voxels ]]; then
        max_slice=$i
        max_num_voxels=$num_voxels
    fi
done
echo $max_slice
cd $cur_dir
rm -rf $tmp_dir
