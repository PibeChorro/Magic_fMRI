#!/usr/bin/bash

# Executable bash file to queue on the headnode to run the freesurfer docker container for reconstruction
# VARIABLES
SUB=${1}
PROJ_DIR=$HOME/Documents/Master_Thesis/DATA/MRI
SUBJECTS_DIR=${PROJ_DIR}/derivatives/freesurfer
COREGISTERED_DIR=${PROJ_DIR}/derivatives/spm12/spm12-preproc/coregistered
meanNiFTI=$COREGISTERED_DIR/$SUB/func/meanu${SUB}_task-magic_bold.nii
# Those are the ROI names you just need to know them -- V1v=1, V1d=2 in the label images (... I guess ???)
# (https://hub.docker.com/r/nben/occipital_atlas) "older version" of neuropythy)

roiname_array=("24dd" "24dv" "a24" "a24pr" "p24" "p24pr" "a32pr" "d32" "p32" "p32pr" "s32"  "33pr" \
"23c" "23d" "d23ab" "31a" "31pd" "31pv" \ 
"44" "45" "47l" "47m" "47s" "a47r" "p47r" \
"AVI" "AAIC" \ 
"IFJa" "IFJp" \
"PHT" "PF") # PHT corresponds to peak values in Danek et al 2014 inferior temporal gyrus and PF to anterior supra marginal gyrus

# ROIs we want to merge
rois_to_merge=("ACC" "PCC" "IFG" "aINSULA" "IFJ")
# all ROIs combined -- this strange syntax is needed

# all ROIs combined -- this strange syntax is needed
allROIs=( "${roiname_array[@]}" "${rois_to_merge[@]}" )

mkdir $SUBJECTS_DIR/$SUB/ROIs

###################################################
# MERGE THE DIFFERENT BRODMAN AREAS TO ONE REGION #
###################################################

# command: mri_mergelabels
# -i an imput image to merge with other input images
# -o output image


# ACC
# LEFT HEMISPHERE
mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/lh.L_24dd_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_24dv_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_a24_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_a24pr_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_p24_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_p24pr_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_a32pr_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_d32_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_p32_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_p32pr_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_s32_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_33pr_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/lh.L_ACC_ROI.label
					
# RIGHT HEMISPHERE
mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/rh.R_24dd_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_24dv_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_a24_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_a24pr_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_p24_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_p24pr_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_a32pr_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_d32_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_p32_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_p32pr_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_s32_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_33pr_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/rh.R_ACC_ROI.label

# PCC 
# LEFT HEMISPHERE 

mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/lh.L_23c_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_23d_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_d23ab_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_31a_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_31pd_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_31pv_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/lh.L_PCC_ROI.label

# RIGHT HEMISPHERE
mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/rh.R_23c_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_23d_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_d23ab_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_31a_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_31pd_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_31pv_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/rh.R_PCC_ROI.label

# IFG 
# LEFT HEMISPHERE
mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/lh.L_44_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_45_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_47l_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_47m_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_47s_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_a47r_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_p47r_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/lh.L_IFG_ROI.label

# RIGHT HEMISPHERE
mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/rh.R_44_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_45_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_47l_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_47m_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_47s_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_a47r_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_p47r_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/rh.R_IFG_ROI.label

# ANTERIOR INSULA
# LEFT HEMISPHERE
mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/lh.L_AVI_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_AAIC_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/lh.L_aINSULA_ROI.label

# RIGHT HEMISPHERE
mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/rh.R_AVI_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_AAIC_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/rh.R_aINSULA_ROI.label

# IFJ
# LEFT HEMISPHERE
mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/lh.L_IFJa_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/lh.L_IFJp_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/lh.L_IFJ_ROI.label

# RIGHT HEMISPHERE
mri_mergelabels -i $SUBJECTS_DIR/${SUB}/label/rh.R_IFJa_ROI.label \
		-i $SUBJECTS_DIR/${SUB}/label/rh.R_IFJp_ROI.label \
		-o $SUBJECTS_DIR/${SUB}/label/rh.R_IFJ_ROI.label


# convert labels into ROIs
for roi in "${allROIs[@]}"
do

	mri_label2vol --label $SUBJECTS_DIR/$SUB/label/lh.L_${roi}_ROI.label \
	--temp $meanNiFTI \
	--reg $SUBJECTS_DIR/$SUB/tkregister_${SUB}.dat \
	--fillthresh 0 \
	--proj frac 0 1 0.1 \
	--subject $SUB \
	--hemi lh \
	--o $SUBJECTS_DIR/$SUB/ROIs/lh.L_Glasser_${roi}.nii

	mri_label2vol --label $SUBJECTS_DIR/$SUB/label/rh.R_${roi}_ROI.label \
	--temp $meanNiFTI \
	--reg $SUBJECTS_DIR/$SUB/tkregister_${SUB}.dat \
	--fillthresh 0 \
	--proj frac 0 1 0.1 \
	--subject $SUB \
	--hemi rh \
	--o $SUBJECTS_DIR/$SUB/ROIs/rh.R_Glasser_${roi}.nii

done
