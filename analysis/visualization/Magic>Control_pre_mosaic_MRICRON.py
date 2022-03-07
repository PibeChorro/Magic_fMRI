###############################################
# CREATE FIGURES OF WHOLE BRAIN UNIVARIATE ANALYSES #
###############################################

# IMPORTANT: This needs to be run in MRIcroGL

# import necessary libraries
# OS interaction and path information
import os
import sys
import glob 
from pathlib import Path
# MRIcroGL functions
import gl 

# get path information
HOME = str(Path.home())
DATA_DIR = os.path.join(HOME, 'Documents', 'Magic_fMRI', 'derivatives', 'spm12', 'spm12-sla')
RESULTS_DIR = os.path.join(HOME, 'ownCloud', 'Magic_scripts_video', 'Paper_Magic_fMRI', 'Images', 'UnivariateResults')
VIDEO_TYPE_DIR = os.path.join(DATA_DIR, 'WholeBrain', 'VideoTypes','9mm-smoothed-mnispace','SpecialMoment')
BACKGROUND_IMG_DIR = os.path.join(HOME,'spm12', 'canonical', 'avg152T1.nii')

# Magic > Control Before revelation
sig_clusters = os.path.join(VIDEO_TYPE_DIR,'Magic > Control Before', '001_30voxthres.nii')
# slices we want to depict (in MNI coordinates)
# axial slices : 60 50 44 30 18 4 -10
# sagital slices = -6
slices_of_interest = "A L+ 60 50 44 30; 18 4 -10 S -6"
# reset defaults
gl.resetdefaults()
# read in background image
gl.loadimage(BACKGROUND_IMG_DIR)
# read in the saved cluster image (manually created in spm)
gl.overlayload(sig_clusters)
gl.colorname(1,'4hot')

gl.mosaic(slices_of_interest )
gl.savebmp(os.path.join(RESULTS_DIR, 'Magic > Control Before001_30voxthres.png'))
#gl.mosaic("A L+ H -0.2 -24 -16 16 40; 48 56 S X R 0");
