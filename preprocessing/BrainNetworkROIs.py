#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# This script shall create binary mask of the seven brain networks worked out 
# by Yeo et al 2011. 

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# The atlas was previously taken from the freesurfer home
# directory and transformed into native space.
# The script iterates over subjects in the freesurfer derivatives directory. 
# Loads in the atlas and creates 7 network ROIs out of it. 
# Each network has a value from 1 - 7 (which value belongs to which network
# was derived from a table showing the color coding of the atlas in the paper)
# The values need to be rounded however because the transformation (probably)
# changed the values slightly.

# interact with the operating system 
import os
import glob
from pathlib import Path
import git
from nilearn.image import new_img_like
import nibabel as nib

################################################
# VARIABLES FOR PATH SELECTION AND DATA ACCESS #
################################################
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
COREGISTERED_DIR= os.path.join(DERIVATIVES_DIR, 'spm12', 'spm12-preproc','coregistered')
FREESURFER_DIR  = os.path.join(DERIVATIVES_DIR,'freesurfer')

SUBJECTS = glob.glob(os.path.join(FREESURFER_DIR,'sub-*'))
SUBJECTS.sort()

NETWORKS = [
    'visual',
    'somatomotor',
    'DAN',
    'VAN',
    'limbic',
    'frontoparietal',
    'DMN']

THRESHOLD = 0.1

for s, sub in enumerate(SUBJECTS):
    # get all permuted accuracy maps of current subject and sort them
    sub_atlas     = glob.glob(os.path.join(FREESURFER_DIR,sub,'atlases','rYeo_7Network.nii'))
    
    img         = nib.load(sub_atlas[0])
    img_data    = img.get_fdata()                  # get data from NIfTI image
    #img_data    = img_data.round()
    img.uncache() 
    
    # iterate over networks and create ROI mask for each
    for n, net in enumerate(NETWORKS):
        output_dir = os.path.join(sub,'corrected_ROIs', NETWORKS[n]+'.nii')
        network_mask = img_data.copy()
        network_mask[(network_mask <= n+1 - THRESHOLD) | (network_mask >= n+1 + THRESHOLD)] = 0
        network_mask[(network_mask > n+1 - THRESHOLD) & (network_mask < n+1 + THRESHOLD)] = 1
        results = new_img_like(ref_niimg=img,data=network_mask)
        nib.save(results, output_dir)
        
        
##################
# WRITE LOG FILE #
##################
# We want to save all important information of the script execution
# To get the git hash we have to check if the script was run locally or on the
# cluster. If it is run on the cluster we want to get the $PBS_O_WORKDIR 
# variable, which preserves the location from which the job was started. 
# If it is run locally we want to get the current working directory.

try:
    script_file_directory = os.environ["PBS_O_WORKDIR"]
except KeyError:
    script_file_directory = os.getcwd()
    
try:
    rep = git.Repo(script_file_directory, search_parent_directories=True)
    git_hash = rep.head.object.hexsha
except git.InvalidGitRepositoryError:
    git_hash = 'not-found'