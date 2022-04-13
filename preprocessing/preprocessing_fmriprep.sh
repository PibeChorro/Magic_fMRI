#!/usr/bin/bash

# Executable bash file to queue on the headnode to run the fmriprep:21.0.1 docker image for preprocessing
SUB=${1}
PROJ_DIR=$HOME/Documents/Magic_fMRI/DATA/MRI
# define the directory where the raw data is stored
CONTAINER_NAME=vplikat_agbartels_fmriprep_${SUB}

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
    --output-spaces func MNI152NLin2009cAsym \
    --write-graph \
    --stop-on-first-crash \
    --fs-no-reconall \
    --no-submm-recon \
    --skip_bids_validation \
    --n_cpus 12 \
    --verbose
