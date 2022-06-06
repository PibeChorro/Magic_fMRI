#!/usr/bin/bash

# Executable bash file to queue on the headnode to run the freesurfer docker container for reconstruction
# VARIABLES
SUB=${1}
PROJ_DIR=$HOME/Documents/Magic_fMRI/DATA/MRI
FS_CONTAINER_NAME=agbartels_freesurfer_${SUB}
SUBJECTS_DIR=/home/derivatives/freesurfer
COREGISTERED_DIR=/home/derivatives/spm12/spm12-preproc_nordic/coregistered
# get the meanEPI image
meanNiFTI=$COREGISTERED_DIR/$SUB/func/meanu${SUB}_task-magic_bold.nii

# RUN THE DOCKER IMAGE
docker run \
--detach \
--name $FS_CONTAINER_NAME \
--rm \
--user "$(id -u):$(id -g)" \
--volume "${PROJ_DIR}:/home" \
--env FS_LICENSE=/home/derivatives/freesurfer/.license \
--env SUBJECTS_DIR=$SUBJECTS_DIR \
--interactive \
--tty \
freesurfer/freesurfer:7.1.1
#docker run: executes docker container create and docker container start
#--detach: container runs in background. If this flag is not set you are stuck in the container terminal and all the following commands are not executed
#--name: A unique name for your container you can use to refer to it (start, stop, restart etc)
#--rm: removes container when stopped
#--user: IMPORTANT - docker runs in root by default. If this flag is not set, every file you create is owned by root and not by yourself
#--volume: kind of mounts a directory into the container
#--env: exports/overwrites an environment variable - here where FS should look for the license file and subjects
#--interactive: this flag is needed so you can send commands to your container
#--tty: opens a terminal in your container

# run docker exec to send the container the actual commands you want to execute
# some commands (like mkdir) can be given to 'docker exec', others however (like cd) need to be given to bash (docker exec bash -c 'bla bla bla').

## create a folder for your subject
#docker exec $FS_CONTAINER_NAME mkdir $SUBJECTS_DIR/${SUB}
## create a subfolder called 'mri'
#docker exec $FS_CONTAINER_NAME mkdir $SUBJECTS_DIR/${SUB}/mri
#create folder for ROI nifti images
docker exec $FS_CONTAINER_NAME mkdir $SUBJECTS_DIR/$SUB/freesurferROIs
## convert your raw anatomical nifti into an mgz file (which needs to be called 001.mgz)
#docker exec $FS_CONTAINER_NAME mri_convert /home/rawdata/${SUB}/anat/${SUB}_T1w.nii /home/derivatives/freesurfer/${SUB}/mri/001.mgz
## execute the recon-all command with the given subject
#docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME recon-all -autorecon-all -subjid ${SUB}

#-------------------------------------------------------------------------------#
# create a registration file
docker exec $FS_CONTAINER_NAME bbregister \
	--mov $meanNiFTI \
	--s $SUB \
	--reg $SUBJECTS_DIR/$SUB/bbregister_${SUB}.dat \
	--bold
# bbregister: registration using a boundary-based cost function
# --mov: template image
# --s: subject id
# --reg: output registration file
# --bold: contrast modality (t1, t2, bold or dti)
#-------------------------------------------------------------------------------#

# create a label and a ROI for the 3rd ventricle as a control ROI
#-------------------------------------------------------------------------------#
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
	mri_vol2label \
	--i $SUB/mri/aseg.mgz \
	--id 15 \
	--l $SUB/label/3rd-ventricle.label

docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
	mri_label2vol \
	--temp $meanNiFTI \
	--label $SUB/label/3rd-ventricle.label \
	--o $SUBJECTS_DIR/$SUB/freesurferROIs/3rd-ventricle.nii \
	--reg $SUBJECTS_DIR/$SUB/bbregister_${SUB}.dat
#-------------------------------------------------------------------------------#

docker stop $FS_CONTAINER_NAME
