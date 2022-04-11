#!/usr/bin/bash

# Executable bash file to queue on the headnode to run the fmriprep:21.0.1 docker image for preprocessing


PROJ_DIR=$HOME/Documents/Magic_fMRI/DATA/MRI
# define the directory where the raw data is stored
RAWDATA=$PROJ_DIR/rawdata
CONTAINER_NAME=vplikat_agbartels_fmriprep_${SUB}

# loop over every subject in RAWDATA and call qsub with docker_fs.sh to reconstruct the subject
for subject in ${RAWDATA}/sub-*
do
    SUB=$(basename "$subject")
    docker run \
        --name $CONTAINER_NAME \
        --interactive \
        --tty \
        --rm \
        --user "$(id -u):$(id -g)" \
        --volume $PROJ_DIR/rawdata:/data:ro \
        --volume $PROJ_DIR/derivatives:/out \
        \
        nipreps/fmriprep:21.0.1 \
        --fs-license-file /out/freesurfer/.license \
        /data /out/fmriprep-21.0.1 \
        participant --participant-label $SUB \
        --output-spaces func \
        --write-graph \
        --stop-on-first-crash \
        --fs-no-reconall \
        --no-submm-recon \
        --skip_bids_validation \
        --verbose
done
