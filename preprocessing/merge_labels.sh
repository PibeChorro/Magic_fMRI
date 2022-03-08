#!/usr/bin/bash

# Executable bash file to queue on the headnode to run the freesurfer docker container to merge the labels of lateral occipital complex LO, ventral occipital complex VO and IPS
# this script is not needed for decoding or other analyses. We create these merged labels for ROI visualization. 
# VARIABLES
SUB=${1}
PROJ_DIR=$HOME/Documents/Magic_fMRI/DATA/MRI
FS_CONTAINER_NAME=agbartels_freesurfer_${SUB}
SUBJECTS_DIR=/home/derivatives/freesurfer

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

#######
# IPS #
#######
# LEFT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
    mri_mergelabels -i ${SUB}/label/lh.wang15atlas.IPS0.label \
			-i ${SUB}/label/lh.wang15atlas.IPS1.label \
			-i ${SUB}/label/lh.wang15atlas.IPS2.label \
			-i ${SUB}/label/lh.wang15atlas.IPS3.label \
			-i ${SUB}/label/lh.wang15atlas.IPS4.label \
			-i ${SUB}/label/lh.wang15atlas.IPS5.label \
			-o ${SUB}/label/lh.wang15atlas.IPS.label 
					
# RIGHT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
	mri_mergelabels -i ${SUB}/label/rh.wang15atlas.IPS0.label \
			-i ${SUB}/label/rh.wang15atlas.IPS1.label \
			-i ${SUB}/label/rh.wang15atlas.IPS2.label \
			-i ${SUB}/label/rh.wang15atlas.IPS3.label \
			-i ${SUB}/label/rh.wang15atlas.IPS4.label \
			-i ${SUB}/label/rh.wang15atlas.IPS5.label \
			-o ${SUB}/label/rh.wang15atlas.IPS.label 
			
######			
# LO #
######
# LEFT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
    mri_mergelabels -i ${SUB}/label/lh.wang15atlas.LO1.label \
                    -i ${SUB}/label/lh.wang15atlas.LO2.label \
                    -o ${SUB}/label/lh.wang15atlas.LO.label
                    
# RIGHT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
    mri_mergelabels -i ${SUB}/label/rh.wang15atlas.LO1.label \
                    -i ${SUB}/label/rh.wang15atlas.LO2.label \
                    -o ${SUB}/label/rh.wang15atlas.LO.label

######
# VO #
######
# LEFT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
    mri_mergelabels -i ${SUB}/label/lh.wang15atlas.VO1.label \
                    -i ${SUB}/label/lh.wang15atlas.VO2.label \
                    -o ${SUB}/label/lh.wang15atlas.VO.label
                    
# RIGHT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
    mri_mergelabels -i ${SUB}/label/rh.wang15atlas.VO1.label \
                    -i ${SUB}/label/rh.wang15atlas.VO2.label \
                    -o ${SUB}/label/rh.wang15atlas.VO.label
                    
docker stop $FS_CONTAINER_NAME
