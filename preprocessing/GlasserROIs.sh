#!/usr/bin/bash

# Executable bash file to queue on the headnode to run the freesurfer docker container for ROI generation
# VARIABLES
SUB=${1}
PROJ_DIR=$HOME/Documents/Magic_fMRI/DATA/MRI
FS_CONTAINER_NAME=agbartels_freesurfer_${SUB}
SUBJECTS_DIR=/home/derivatives/freesurfer
COREGISTERED_DIR=/home/derivatives/spm12/spm12-preproc/coregistered
# get the meanEPI image
meanNiFTI=$COREGISTERED_DIR/$SUB/func/meanu${SUB}_task-magic_bold.nii
bilateral_rois=("IFJ" "44" "6r" "BA6" "FEF" "pACC" "mACC" "aACC" "8BM" "AI" "AVI" "PH")
left_lateral_rois=("IFS" "45" "BA46" "BA8" "BA9" "IPC")

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

mkdir $SUBJECTS_DIR/$SUB/ROIs

###########################################
# MERGE THE DIFFERENT LABELS IF NECESSARY #
###########################################

# command: mri_mergelabels
# -i an imput image to merge with other input images
# -o output image

# BILATERAL ROIS

# INFERIOR FRONTAL GYRUS BILATERAL FROM DANEK ET AL 2015 WHOLE VIDEO
# x   y   z
# -44 4   18
# 44  12  24
# joined in IFJ (IFJa, IFJp), BA 44 and BA 6r

# only IFJ has to be merged
# LEFT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_IFJa_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_IFJp_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_IFJ_ROI.label

# RIGHT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_IFJa_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_IFJp_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/rh.R_IFJ_ROI.label

# SUPERIOR/MIDDLE FRONTAL GYRUS BILATERAL FROM DANEK ET AL 2015 MAGIC MOMENT
# x   y   z
# -24	8	52
# 28	8	52
# joined in BA6 (6a, 6ma, i6-8) and FEF

# only BA6 has to be merged
# LEFT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_6a_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_6ma_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_i6-8_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_BA6_ROI.label
# RIGHT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_6a_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_6ma_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_i6-8_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/rh.R_BA6_ROI.label

# ANTERIOR CINGULATE CORTEX BILATERAL FROM DANEK ET AL 2015 MAGIC MOMENT and PARRIS ET AL 2009
# x   y   z
# -4	28	40 (Danek)
# 0	  0	  26 (Danek)
# -4	38	19 (Parris)
# joined in pACC (p24pr a24pr 33pr p32pr), mACC (d32 a32pr p24), aACC (a24 p32 s32) and BA 8BM

# only pACC, mACC and aACC have to be merged
# LEFT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_p24pr_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_a24pr_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_33pr_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_p32pr_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_pACC_ROI.label

docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_d32_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_a32pr_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_p24_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_mACC_ROI.label

docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_a24_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_p32_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_s32_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_aACC_ROI.label

# RIGHT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_p24pr_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_a24pr_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_33pr_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_p32pr_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/rh.R_pACC_ROI.label

docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_d32_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_a32pr_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_p24_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/rh.R_mACC_ROI.label

docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_a24_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_p32_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_s32_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/rh.R_aACC_ROI.label


# ANTERIOR INSULA BILATERAL FROM DANEK ET AL 2015 - MAGIC MOMENT
# x   y   z
# -32	20	-4
# 34	20	-2
# joined by AVI and AI (AAIC, MI)
# LEFT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_MI_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_AAIC_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_AI_ROI.label

# RIGHT HEMISPHERE
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_MI_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/rh.R_AAIC_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/rh.R_AI_ROI.label

# INFERIOR TEMPORAL GYRUS, TEMPORO-OCCIPITAL DIVISION BILATERAL BY DANEK ET AL 2015 - MAGIC MOMENT
# x   y   z
# 62	-56	-8
# -46	-60	-12
# joined by PH -> no merging necessary


#####################
# LEFT LATERAL ROIS #
#####################

# INFERIOR FRONTAL GYRUS PARS TRINGULARIS LEFT FROM DANEK ET AL 2015 MAGIC MOMENT
# x   y   z
# -52	34	10
# joined in IFS (IFSa, IFSp), BA 45, BA 46 (a9-46v, p9-46v 46)

# only IFS and BA 46 have to be merged
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_IFSa_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_IFSp_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_IFS_ROI.label

docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_a9-46v_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_p9-46v_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_46_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_BA46_ROI.label

# MIDDLE FRONTAL GYRUS LEFT FROM PARRIS ET AL 2009
# x   y   z
# -22	36	44
# joined in BA8 (8C 8Ad 8BL 8Av)

docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_8C_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_8Ad_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_8BL_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_8Av_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_BA8_ROI.label


# LEFT BA9 ???
# joined in BA9 (9ü, 9-46d)
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_9p_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_9-46d_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_BA9_ROI.label

# ANTERIOR SUPRAMARGINAL GYRUS LEFT BY DANEK ET AL 2015 - MAGIC MOMENT
# x   y   z
# -66	-32	32
# joined in inferior parietal cortex - IPC (PF, PFt, PFop)
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
  mri_mergelabels \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_PF_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_PFt_ROI.label \
    -i $SUBJECTS_DIR/${SUB}/label/lh.L_PFop_ROI.label \
    -o $SUBJECTS_DIR/${SUB}/label/lh.L_IPC_ROI.label

# command: mri_label2vol
# --label input label image, that you want to convert
#	--temp template volume - the output has the same size and geometry
#	--reg a registration matrix file
#	--fillthresh a value that needs to be exceeded in order for a voxel to be considered a member of the label
#	--proj Type can be either abs or frac. abs means that the start, stop, and delta are measured in mm. frac means that start, stop, and delta are relative to the thickness at each vertex. The label definition is changed to fill in label points in increments of delta from start to stop. Uses the white surface. The label MUST have been defined on the surface.
#	--subject the subject in your freesurfer subject directory
#	--hemi which hemisphere (lh or rh)
#	--o output file

# first convert bilateral labels into ROIs
for roi in "${bilateral_rois[@]}"
do

  docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
    mri_label2vol --label $SUBJECTS_DIR/$SUB/label/lh.L_${roi}_ROI.label \
    --temp $meanNiFTI \
    --reg $SUBJECTS_DIR/$SUB/tkregister_${SUB}.dat \
    --fillthresh 0 \
    --proj frac 0 1 0.1 \
    --subject $SUB \
    --hemi lh \
    --o $SUBJECTS_DIR/$SUB/ROIs/lh.L_Glasser_${roi}.nii
  docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
    mri_label2vol --label $SUBJECTS_DIR/$SUB/label/rh.R_${roi}_ROI.label \
    --temp $meanNiFTI \
    --reg $SUBJECTS_DIR/$SUB/tkregister_${SUB}.dat \
    --fillthresh 0 \
    --proj frac 0 1 0.1 \
    --subject $SUB \
    --hemi rh \
    --o $SUBJECTS_DIR/$SUB/ROIs/rh.R_Glasser_${roi}.nii
done

# second only convert left lateral labels
for roi in "${left_lateral_rois[@]}"
do
  docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME \
    mri_label2vol --label $SUBJECTS_DIR/$SUB/label/lh.L_${roi}_ROI.label \
    --temp $meanNiFTI \
    --reg $SUBJECTS_DIR/$SUB/tkregister_${SUB}.dat \
    --fillthresh 0 \
    --proj frac 0 1 0.1 \
    --subject $SUB \
    --hemi lh \
    --o $SUBJECTS_DIR/$SUB/ROIs/lh.L_Glasser_${roi}.nii
done

docker stop $FS_CONTAINER_NAME