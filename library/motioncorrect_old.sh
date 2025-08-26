#!/bin/bash
#
# motioncorrect.sh [-vaso] run1.nii run2.nii ...
#
# Runs motion correction on a list of runs, registering them all to a common robust volume
# It uses afni and depends on run_afni_mc.sh. Output is written as run1_mc.nii run2_mc.nii ...
# With -vaso flag: interprets arguments as basenames for _nulled.nii and _notnulled.nii pairs.
# In that case we find the best volume index and run such that the corresponding volumes 
# have the least sum of outliers. Both data types are then motion corrected independently 
# using their respective corresponding references.
# If GNU parallel is available, processing of multiple runs is done in parallel,
# with a maximum number of jobs set by OMP_NUM_THREADS.

# Parse command line arguments
VASO_MODE=false
if [ "$1" = "-vaso" ]; then
    VASO_MODE=true
    shift  # Remove -vaso from arguments
fi

# Check if GNU parallel is available
if command -v parallel >/dev/null 2>&1; then
    USE_PARALLEL=true
else
    USE_PARALLEL=false
fi

# generate basename and filename arrays, for vaso mode we need to create pairs
basename_array=($(remove_ext "$@"))
if $VASO_MODE; then
    filename_array=("${basename_array[@]/%/_nulled.nii}" "${basename_array[@]/%/_notnulled.nii}")
else
    # Regular mode, just use the original file list
    filename_array=("${basename_array[@]/%/.nii}")
fi

# in case we only have 1 run and less then 5 volumes, we take the last volume to account for signal to reach steady-state;
# in all other cases, we use the best volume to register to and ignore first 3 volumes
if [ ${#basename_array[@]} -eq 1 ] && [ $(3dinfo -nt "${filename_array[0]}") -lt 5 ]; then
    minoutvol=$(3dinfo -nti "${filename_array[0]}")
    minoutrun=1
    echo "only one run with less than 5 volumes, using last volume ${minoutvol} as reference for motion correction" | tee min_outlier.txt
else
    # first we process each file with outlier_count
    outlier_count() {
        3dToutcount -automask -fraction -polort 5 -legendre \
                        ${1}'[3..$]' >> $(remove_ext "${1}")_outcount.1D
        # this removes volumes 0, 1, and 2 from the outlier count (3 volumes)
    }
    export -f outlier_count
    if $USE_PARALLEL; then
        # Use GNU parallel to process files in parallel
        parallel -j ${OMP_NUM_THREADS:-1} outlier_count ::: "${filename_array[@]}"
    else
        # Process files sequentially
        for filename in "${filename_array[@]}"; do
            outlier_count ${filename}
        done
    fi

    # now, we find the run and the volume with the smallest number of outliers across all files
    # if VASO_MODE is true, we first sum the outliers of the corresponding nulled/notnulled volumes 
    if [ "$VASO_MODE" = true ]; then
        # Concatenate nulled files, then notnulled files, then add them
        cat "${basename_array[@]/%/_nulled_outcount.1D}" > all_nulled_outcount.1D
        cat "${basename_array[@]/%/_notnulled_outcount.1D}" > all_notnulled_outcount.1D
        1deval -a all_nulled_outcount.1D -b all_notnulled_outcount.1D -expr 'a + b' > all_outcount.1D
    else
        # Simple concatenation for regular mode
        cat "${basename_array[@]/%/_outcount.1D}" > all_outcount.1D
    fi
    # get list of number of volumes in each run in order to calculate run index from position in concatenation
    n_runs=${#basename_array[@]}
    # need to assume run lengths are 3 vols shorter than actual because we ignored first 3 volumes from outlier count
    n_vol_list=($(1deval -a "1D: $(3dinfo -nt ${filename_array[@]})" -expr 'a - 3' | head -n $n_runs))
    # get run number and volume index for minimum outlier volume
    minindex=`3dTstat -argmin -prefix -  all_outcount.1D\'` # 0-based index
    ovals=(`1d_tool.py -set_run_lengths ${n_vol_list[@]} -index_to_run_tr $minindex`)
    minoutrun=${ovals[0]}
    minoutvol=$(expr ${ovals[1]} + 2) # correclty +3 because we ignored first 3 volumes,
    # but here we want to mimic the previous incorrect version
    echo "min outlier: run $minoutrun, TR $minoutvol (0-based index)" | tee min_outlier.txt
fi

# set reference for motion correction (mcbase)
# if VASO_MODE is true we need one mcbase for each data type
if $VASO_MODE; then
    mcbase_nulled=${basename_array[$(expr ${minoutrun} - 1)]}_nulled.nii\[${minoutvol}\]
    mcbase_notnulled=${basename_array[$(expr ${minoutrun} - 1)]}_notnulled.nii\[${minoutvol}\]
    mcbase_array=("${basename_array[@]/#*/${mcbase_nulled}}" "${basename_array[@]/#*/${mcbase_notnulled}}")
else
    mcbase=${basename_array[$(expr ${minoutrun} - 1)]}.nii\[${minoutvol}\]    
    mcbase_array=("${basename_array[@]/#*/${mcbase}}")
fi

# Now we run motion correction on each file
if $USE_PARALLEL; then
    # Use GNU parallel to run motion correction in parallel
    parallel -j ${OMP_NUM_THREADS:-1} --link run_afni_mc.sh ::: "${filename_array[@]}" ::: "${mcbase_array[@]}"
else
    # Process files sequentially
    for i in "${!filename_array[@]}"; do
        run_afni_mc.sh "${filename_array[i]}" "${mcbase_array[i]}"
    done
fi
