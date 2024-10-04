#!/usr/bin/env bash

# Generates a FreeSurfer group average subject from scratch (without fsaverage)
# Denis Chaimow 2024-01-24

set -Eeuo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)

usage() {
    cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") -d <SUBJECTS_DIR> [Options]

Required:
    -d <SUBJECTS_DIR>                  Dirrectory containing processed FreeSurfer data sets.
                                        All subdirectorie matching 'sub-\*' will be used.
Options:
    -s <talairach|group|TEMPLATE_FILE> Specify which common volume space to use. 
                                        Can be either 'talairach' to use talairach.xfm from 
                                        recon-all, or 'group' to calculate a within group 
                                        registration space, or a template file can be provided 
                                        to which all subjects will be registered.
                                        (default: talairach)
    -n <AVGSUBNAME>                    Name of average subject to be created.
                                        (default: myaverage)
    -p n                               Use GNU parallel (if available) with n parallel jobs
    -v                                 Verbose
    -h                                 Display help
EOF
}

cleanup() {
  trap - SIGINT SIGTERM ERR EXIT
  # script cleanup here
  cd ${curDir}
  rm -rf ${tmpdir}
}

# default parameters
SPACE=talairach
AVGSUBNAME=myaverage
PARALLEL=false
SUBJECTS_DIR=
NJOBS=

# parse command line parameters
while getopts ":d:s:n:p:vh" opt; do
    case ${opt} in
        d) SUBJECTS_DIR=${OPTARG};;
        s) SPACE=${OPTARG};;
        n) AVGSUBNAME=${OPTARG};;
        p) if command -v parallel; then NJOBS=${OPTARG}; fi;;
        v) set -x;;
        h) usage; exit ;;
        :) usage 1>&2; exit;;
        *) usage 1>&2; exit;;
    esac
done
if [ -z ${SUBJECTS_DIR} ]; then usage 1>&2; exit; fi

# generate list (array) of subjects
SUBJECTS=( $(find $SUBJECTS_DIR -iname "sub-*" -type d | xargs -L1 -I{} basename "{}") )


# switch to temporary directory
curDir=$(pwd)
tmpdir=$(mktemp -d)
cd ${tmpdir}
trap cleanup SIGINT SIGTERM ERR EXIT

# check which volume space to use and perform registrations if necessary
if [[ -f $SPACE ]]; then
    # if existing template file provided, register to that
    TEMPLATE_FILE=$SPACE
    SPACE=$(basename $(basename $(basename $SPACE .mgz) .nii) .nii.gz)
    if [ -n ${NJOBS} ]; then
        parallel --jobs $NJOBS \
                 mri_coreg \
                 --mov $SUBJECTS_DIR/{}/mri/norm.mgz \
                 --ref $TEMPLATE_FILE \
                 --reg $SUBJECTS_DIR/{}/mri/transforms/${SPACE}.lta \
                 --dof 12 ::: ${SUBJECTS[@]}
    else
        for SUBJECT in ${SUBJECTS[@]}; do
            mri_coreg \
                 --mov $SUBJECTS_DIR/${SUBJECT}/mri/norm.mgz \
                 --ref $TEMPLATE_FILE \
                 --reg $SUBJECTS_DIR/${SUBJECT}/mri/transforms/${SPACE}.lta \
                 --dof 12 
        done
    fi                   
elif [[ $SPACE == "group" ]]; then
    # if group volume space is required for averaging
    # then first construct an unbiased robust template and save registrations
    SUBJECTS_SUBDIRS=("${SUBJECTS[@]/#/${SUBJECTS_DIR}\/}")
    mri_robust_template --mov ${SUBJECTS_SUBDIRS[@]/%/\/mri/norm.mgz} \
                        --lta ${SUBJECTS_SUBDIRS[@]/%/\/mri/transforms/group.lta} \
                        --template group_template.mgz \
                        --inittp 0 \
                        --affine \
                        --satit
elif [[ $SPACE != "talairach" ]]; then
     echo "SPACE must be one of talairach, group or a template file!" >1&2
     exit
fi
# otherwise simply use existing talairach space registrations from recon-all

for hemi in lh rh; do
    # create "seed" template from sub-01
    mris_make_template ${hemi} sphere ${SUBJECTS[0]} ${hemi}.mytemp0.tif
    # results in single template file ?h.mytemp0.tif 
    
    # register all subjects to mytemp0.tif
    if [ -n ${NJOBS} ]; then
        parallel --jobs $NJOBS \
                 mris_register $SUBJECTS_DIR/{}/surf/${hemi}.sphere ${hemi}.mytemp0.tif \
                 $SUBJECTS_DIR/{}/surf/${hemi}.sphere.myreg0 ::: ${SUBJECTS[@]}
    else
        for SUBJECT in ${SUBJECTS[@]}; do
             mris_register $SUBJECTS_DIR/${SUBJECT}/surf/${hemi}.sphere ${hemi}.mytemp0.tif \
                           $SUBJECTS_DIR/${SUBJECT}/surf/${hemi}.sphere.myreg0
        done
    fi
    # results in transformed spherical surface for each subjects ?h.sphere.myreg0
    
    # next round: create new template from all already registered subjects
    mris_make_template ${hemi} sphere.myreg0 ${SUBJECTS[@]} ${hemi}.mytemp1.tif
    # results in single template file ?h.mytemp1.tif 
    
    
    # register all subjects to mytemp1.tif
    if [ -n ${NJOBS} ]; then
        parallel --jobs $NJOBS \
	         mris_register $SUBJECTS_DIR/{}/surf/${hemi}.sphere ${hemi}.mytemp1.tif \
                 $SUBJECTS_DIR/{}/surf/${hemi}.sphere.myreg ::: ${SUBJECTS[@]}
    else
        for SUBJECT in ${SUBJECTS[@]}; do
             mris_register $SUBJECTS_DIR/${SUBJECT}/surf/${hemi}.sphere ${hemi}.mytemp1.tif \
                           $SUBJECTS_DIR/${SUBJECT}/surf/${hemi}.sphere.myreg
        done
    fi
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

if [[ $SPACE==talairach ]]; then
    make_average_subject --out $AVGSUBNAME \
                         --surf-reg sphere.myreg \
                         --sd ${SUBJECTS_DIR}\
                         --subjects ${SUBJECTS[@]}
else
    make_average_subject --out $AVGSUBNAME \
                         --surf-reg sphere.myreg \
                         --xform ${SPACE}.lta \
                         --sd ${SUBJECTS_DIR} \
                         --subjects ${SUBJECTS[@]} 
fi

cp lh.mytemp.tif ${curDir}/lh.${AVGSUBNAME}.tif
cp rh.mytemp.tif ${curDir}/rh.${AVGSUBNAME}.tif
