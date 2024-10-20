#!/bin/bash
# runs command on multiple subjects (in parallel if possible)

usage () {
    echo "usage: $(basename $0) -c command -d studyDataDir [-l subjList.txt] [-n nJobs]"
    echo "  -c command: command to run on each subject"
    echo "  -d studyDataDir: directory containing subject data"
    echo "  -l subjList.txt: file containing list of subjects to process"
    echo "  -n nJobs: number of jobs to run in parallel"
}

# 1. parse command line arguments
while getopts "c:d:l:n:h" opt; do
    case $opt in
        c)
            command=$OPTARG;;
        d)
            studyDataDir=$OPTARG;;
        l)
            subjList_file=$OPTARG;;
        n)
            nJobs=$OPTARG;;
        h)
            usage;;
        *)
            usage
            exit 1;;
    esac
done

# 2. check if all required arguments are provided
if [ -z "$command" ] || [ -z "$studyDataDir" ]
then
    usage
    exit 1
fi

# 3. generate subject list
if [ ! -f "$subjList_file" ]
then
    subjectList=$(ls ${studyDataDir} | grep ^sub-)
else
    subjectList=$(cat $subjList_file)
fi

# 4. check if GNU parallel is installed
if [ -z $(which parallel) ]
then
    echo "GNU parallel is not installed, running the command sequentially"
    for subject in $subjectList
    do
        $command $studyDataDir $subject
    done
    exit 0
fi  

# 5. determine the number of jobs to run in parallel
nSubjects=$(echo $subjectList | wc -w)
nCores=$(nproc)

# if njobs is not provided, run one job per subject but not more than the number of cores
if [ -z "$nJobs" ]
then
    nJobs=$(($nSubjects<$nCores?$nSubjects:$nCores))
fi 
echo "Running $nSubjects subjects using $nJobs jobs in parallel"

# 6. run the command in parallel
parallel --jobs $nJobs $command $studyDataDir {} ::: $subjectList
