#!/bin/bash
# usage: ./create_fs_group_average_subject.sh <SUBJECTS_DIR> <group|talairach> <avgerage_subject_name>
# generates a fs group average subject from scratch

SUBJECTS_DIR=$1
SPACE=$(tr '[:upper:]' '[:lower:]' <<< $2)
AVGSUBNAME=${3:-myaverage}
NJOBS=8

SUBJECTS=( $(find $SUBJECTS_DIR -iname "sub-*" -type d | xargs -L1 -I{} basename "{}") )

curDir=$(pwd)
tmpdir=$(mktemp -d)
cd ${tmpdir}

# if group volume space is required for averaging
# then first construct an unbiased robust template and save registrations
if [[ $SPACE == "group" ]]; then
    SUBJECTS_SUBDIRS=("${SUBJECTS[@]/#/${SUBJECTS_DIR}\/}")
    mri_robust_template --mov ${SUBJECTS_SUBDIRS[@]/%/\/mri/norm.mgz} \
                        --lta ${SUBJECTS_SUBDIRS[@]/%/\/mri/transforms/group.lta} \
                        --template group_template.mgz \
                        --satit \
                        --affine
fi

for hemi in lh rh; do
    # create "seed" template from sub-01
    mris_make_template ${hemi} sphere ${SUBJECTS[0]} ${hemi}.mytemp0.tif
    # results in single template file ?h.mytemp0.tif 
    
    # register all subjects to mytemp0.tif 
    parallel --jobs $NJOBS \
             mris_register $SUBJECTS_DIR/{}/surf/${hemi}.sphere ${hemi}.mytemp0.tif \
             $SUBJECTS_DIR/{}/surf/${hemi}.sphere.myreg0 \
             ::: ${SUBJECTS[@]}
    # results in transformed spherical surface for each subjects ?h.sphere.myreg0
    
    # next round: create new template from all already registered subjects
    mris_make_template ${hemi} sphere.myreg0 ${SUBJECTS[@]} ${hemi}.mytemp1.tif
    # results in single template file ?h.mytemp1.tif 
    
    
    # register all subjects to mytemp1.tif
    parallel --jobs $NJOBS \ 
	     mris_register $SUBJECTS_DIR/{}/surf/${hemi}.sphere ${hemi}.mytemp1.tif \
             $SUBJECTS_DIR/{}/surf/${hemi}.sphere.myreg \
             :::  ${SUBJECTS[@]}
    # results in transformed spherical surface for each subjects ?h.sphere.myreg1
    
    # finally: create final template from all already registered subjects
    mris_make_template ${hemi} sphere.myreg ${SUBJECTS[@]} ${hemi}.mytemp.tif
    # results in single template file ?h.mytemp.tif 
    
    # repeat everything for the other hemisphere:
done

# the result of the above processing is twofold:
# 1. a template, which is a single .tif file, which contains information about average and
#    variance of sulcal patterns as a function of spherical coordinates
#    - inflated.h = mean curvature of inflated surface
#    - curv = spatially smoothed mean curvature of white surface
#    - sulc = spatially smoothed convecity of white surface
# 2. for each subject a spherical surface such that for every vertex of the "original" surfaces
#    it gives us the spherical coordinates of the corresponding vertex providing a link to the
#    common spherical coordinate system

# We can then use 1 (the template file) to register (using mris_register) other subjects to the
# same spherical coordinate system based on correspondence of sulcal patterns.
#
# From 2 (spherical transformed or registered surface) and all the other info in the fs
# directories of individual subjects we can then generate an average subject.

make_average_subject --out $AVGSUBNAME --surf-reg sphere.myreg --xform ${SPACE}.lta --subjects ${SUBJECTS[@]} 

cp lh.mytemp.tif ${curDir}/lh.${AVGSUBNAME}.tif
cp rh.mytemp.tif ${curDir}/rh.${AVGSUBNAME}.tif

cd ${curDir}
