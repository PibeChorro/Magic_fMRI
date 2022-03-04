#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########################
# PURPOSE OF THIS SCRIPT #
##########################
# The purpose of this script is to create a group mask, so that in further 
# analyses there are no subject specific voxels 

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# FIRST STEP 
# Get all important path information and information about the analysis from 
# command line input.
# SECOND STEP 
# Get all subjects in the GLM directory asked and read in the mask of every
# subject. Add them all up to one large np.matrix and turn it into a bool map
# by comparing it to the number of subjects
# THIRD STEP
# Save group mask

#############
# LIBRARIES #
#############

# FIRST STEP
# interact with the operating system 
import os
import glob
import argparse
from pathlib import Path
import git
# data structuration and calculations
import numpy as np   # most important numerical calculations
# read in mat files
from nilearn.image import new_img_like
import nibabel as nib
# optimize time performance
import time

# get start time
T_START = time.time()

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--smooth",             nargs='?',  const=0,        
                    default=0,      type=int)   # what data should be used
parser.add_argument("--analyzed",           nargs='?', const='moment',  
                    default='video',   type=str)
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
SMOOTHING_SIZE  = ARGS.smooth
ANALYZED        = ARGS.analyzed

if ANALYZED == 'moment':
    data_analyzed = 'SpecialMoment'
elif ANALYZED == 'video':
    data_analyzed = 'WholeVideo'
else:
    raise
    
################################################
# VARIABLES FOR PATH SELECTION AND DATA ACCESS #
################################################
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
# where to look for the beta images
if SMOOTHING_SIZE > 0:
    GLM_DATA_DIR    = str(SMOOTHING_SIZE)+'mm-smoothed-mnispace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               data_analyzed)
else:
    GLM_DATA_DIR    = 'mnispace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               data_analyzed)
    
OUTPUT_DIR = os.path.join(FLA_DIR, 'group_mask.nii')

# SECOND STEP
SUBJECTS = glob.glob(os.path.join(FLA_DIR,'sub-*'))
SUBJECTS.sort()

masks = []

for s, sub in enumerate(SUBJECTS):
    # get all permuted accuracy maps of current subject and sort them
    sub_mask     = glob.glob(os.path.join(FLA_DIR,sub,'mask.nii'))
    
    # read in the randomly selected accuracy map, transform into numpy array and add it to the list
    img             = nib.load(sub_mask[0])
    img_data        = img.get_fdata()                  # get data from NIfTI image
    img.uncache() 
    masks.append(img_data)
    
masks = np.array(masks)
group_mask = masks.sum(axis=0)
group_mask = group_mask==len(SUBJECTS)

# THIRD STEP
results = new_img_like(ref_niimg=img,data=group_mask)
nib.save(results, OUTPUT_DIR)