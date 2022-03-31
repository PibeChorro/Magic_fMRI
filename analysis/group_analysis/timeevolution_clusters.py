#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# The data for this script is derived from an fMRI experiment in which subjects
# viewed different videos. Either magic, control or surprise videos.
# The experiment was devided into 3 blocks. Each block consisted of 4
# experimental runs. In each run subjects viewed 24 videos (each video is
# considered a trial).
# The videos in each block were associated with one object (Balls, Cards and
# Sticks) and there were 3 magic effects (Appear, Change and Vanish). For each
# magic effect and object there are two trick versions (i.e. Appear1, Appear2,
# Change1,...). This resulted in 6 magic videos per object = 18 magic videos
# and for every magic video there was a corresponding control video showing
# the same movements without the magical effect. Additionally per object there
# were 3 surprise videos showing unusual surprising actions performed with the
# objects (e.g. eating a playing card).
# After the second run in each block the underlying method behind each magic
# trick was presented.

#                       TIME
#   ---------------------------------->
#   OBJECT1 R1  R2  Revelation  R3  R4  |
#   OBJECT2 R1  R2  Revelation  R3  R4  |   TIME
#   OBJECT3 R1  R2  Revelation  R3  R4  v

# RUNS: 2*Appear1 Magic 2*Appear2 Magic 2*Appear1 Control 2*Appear2 Control
#       2*Vanish1 Magic 2*Vanish2 Magic 2*Vanish1 Control 2*Vanish2 Control
#       2*Change1 Magic 2*Change2 Magic 2*Change1 Control 2*Change2 Control
#       2*Surprise1     2*Surprise2     2*Surprise13
#       = 24 Videos

# The aim of the experiment was to find neural correlates of surprise and in
# particular of surprising events we consider "impossible".
# The data used are beta estimate NIfTI images derived from a GLM using SPM12
# in MATLAB.

##########################
# PURPOSE OF THIS SCRIPT #
##########################
# The purpuse of this script is to create ROIs from significant clusters in other analyses.
# Specifically for significant clusters in the univariate whole brain non-parametric contrasts
# The resulting nifti images from SnPM are images containing the negative log with basis 10 of the p-value:
# -log_10(0.1)=1    -log_10(0.01)=2    -log_10(0.001)=3
# Connected voxels from these maps with a value >= 3 are considered to be a ROI.

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# Get all important libraries, and paths
# Then read in the p-value maps in and identify all connected voxels, that have a value equal or higher than 3
# Check the size of the cluster. If it is large enough, save it as a ROI

#############
# LIBRARIES #
#############

# interact with the operating system
import glob
import os
import argparse
import warnings
from pathlib import Path
# data structuration and calculations
import pandas as pd
from nilearn.image import load_img, new_img_like
import nibabel as nib
import numpy as np  # most important numerical calculations
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
# library for neuroimaging
# optimize time performance
import time

# make convergense warning an error in oder to catch it
warnings.simplefilter('error', ConvergenceWarning)

# get start time
T_START = time.time()

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--smooth", nargs='?', const=0,
                    default=0, type=int)  # what data should be used
parser.add_argument("--time", nargs='?', const='moment',
                    default='moment', type=str)  # what data should be used
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values
TIME = ARGS.time

if TIME == 'moment':
    time_analyzed = 'SpecialMoment'
elif TIME == 'video':
    time_analyzed = 'WholeVideo'

# variables for path selection and data access
HOME = str(Path.home())
PROJ_DIR = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
SPM_FLA_DIR = os.path.join(DERIVATIVES_DIR, 'spm12', 'spm12-fla', 'WholeBrain')
TYPE_RES_DIR = os.path.join(DERIVATIVES_DIR, 'snpm13', 'snpm13-sla', 'WholeBrain', 'VideoTypes',
                            '6mm-smoothed-mnispace', time_analyzed)
EFFECT_RES_DIR = os.path.join(DERIVATIVES_DIR, 'snpm13', 'snpm13-sla', 'WholeBrain', 'MagicEffects',
                              '6mm-smoothed-mnispace', time_analyzed)

ALL_RUNS = np.linspace(1, 12, 12, dtype=int)

TYPE_CONTRASTS_OF_INTEREST = ['Magic Before vs Magic After', 'MagPre-ConPre vs MagPost-ConPost']
EFFECT_CONTRASTS_OF_INTEREST = ['Appear Before vs Appear After', 'AppPre-ConPre vs AppPost-ConPost',
                                'Change Before vs Change After', 'ChaPre-ConPre vs ChaPost-ConPost',
                                'Vanish Before vs Vanish After', 'VanPre-ConPre vs VanPost-ConPost']
EFFECTS = ['Appear', 'Appear',
           'Change', 'Change',
           'Vanish', 'Vanish']

CLUSTER_IMAGE = 'clusters.nii'
CLUSTERSIZE_THRESHOLD = 10

for type_contrast in TYPE_CONTRASTS_OF_INTEREST:

    # dictionary saving the AIC values for each model
    aic_dict = {
        'cluster': [],
        'pre_post': [],
        'runs_in_block': [],
        'all_runs': []
    }

    # FIRST DO MAGIC BEFORE VS MAGIC AFTER
    img_path = os.path.join(TYPE_RES_DIR, type_contrast, CLUSTER_IMAGE)
    # read in a dataframe containing the video condition, the prepost condition and the corresponding beta names
    # created by 'create_type_ROI_DF.py'
    df_path = os.path.join(DERIVATIVES_DIR, 'univariate-ROI', 'SpecialMoment', 'VideoType', 'data_frame.csv')
    dataframe = pd.read_csv(df_path)
    # reduce dataframe to only Magic conditions
    magic_df = dataframe[dataframe.Type == 'Magic']

    # get all subjects
    subjects_path = os.path.join(SPM_FLA_DIR, 'VideoTypes', '6mm-smoothed-mnispace', time_analyzed)
    subjects = glob.glob(os.path.join(subjects_path, 'sub-*'))
    subjects.sort()
    ALL_SUBJECTS = np.linspace(1, len(subjects), len(subjects), dtype=int)

    # create variables that are the same for all clusters
    subs = np.repeat(ALL_SUBJECTS, len(ALL_RUNS))
    runs = np.tile(ALL_RUNS, len(subjects))
    pre_post_inblock    = (runs-1)//2   # Resulting in 0,0,1,1,2,2,3,3 ...
    # pre_post_inblock%2 resulting in 0,0,1,1,0,0,1,1 ...
    pre_post = ['pre' if pp%2 == 0 else 'post' for pp in pre_post_inblock]
    runs_in_block = (runs-1) % 4        # Resulting in 0,1,2,3,0,1,2,3, ...

    current_img = load_img(img_path)
    img_data = current_img.get_fdata()

    # get the highest value in our clusters image - to get the number of separate clusters
    num_clusters = np.max(img_data, axis=None).astype(int)

    for c in range(num_clusters):
        # set start parameter for model. None if it is the first time, else take the values from the previous model
        if c == 0:
            params_pre_post = None
            params_runs_in_block = None
            params_all_runs = None
        else:
            params_pre_post = pre_post_res.params
            params_runs_in_block = runs_in_block_res.params
            params_all_runs = runs_res.params
        # get cluster mask
        cluster_mask = np.array(img_data == c + 1)
        # check wether the cluster is large enough
        if np.sum(cluster_mask, axis=None) < CLUSTERSIZE_THRESHOLD:
            continue

        # fill dictionary with entries
        data_dict = {
            'subject': subs,
            'runs': runs,
            'pre_post': pre_post,
            'runs_in_block': runs_in_block,
            'betas': []
        }

        for s, sub in enumerate(subjects):
            sub_df = magic_df[magic_df.subject == s+1]
            for current_run in ALL_RUNS:
                current_beta = sub_df.BetaNames[sub_df.Runs == current_run].values[0]
                beta_nii = nib.load(os.path.join(sub, current_beta))
                beta_data = beta_nii.get_fdata()
                beta_roi = beta_data[cluster_mask & ~np.isnan(beta_data)]

                data_dict['betas'].append(beta_roi.mean(axis=None))

        cluster_df = pd.DataFrame(data=data_dict, columns=data_dict.keys())
        cluster_df.to_csv(path_or_buf=os.path.join(TYPE_RES_DIR, type_contrast, 'cluster{:02}_df.csv'.format(c+1)), index=False)

        aic_dict['cluster'].append(c+1)

        pre_post_model = smf.mixedlm("betas ~ pre_post", cluster_df, groups=cluster_df['subject'])
        try:
            pre_post_res = pre_post_model.fit(reml=False, method="cg")
        except ConvergenceWarning:
            pre_post_res = pre_post_model.fit(reml=False, method="cg", start_params=pre_post_res.params)
        aic_dict['pre_post'].append(pre_post_res.aic)

        runs_in_block_model = smf.mixedlm("betas ~ runs_in_block", cluster_df, groups=cluster_df['subject'])
        try:
            runs_in_block_res = runs_in_block_model.fit(reml=False, method="cg")
        except ConvergenceWarning:
            runs_in_block_res = runs_in_block_model.fit(reml=False, method="cg", start_params=runs_in_block_res.params)
        aic_dict['runs_in_block'].append(runs_in_block_res.aic)

        runs_model = smf.mixedlm("betas ~ runs", cluster_df, groups=cluster_df['subject'])
        try:
            runs_res = runs_model.fit(reml=False, method="cg")
        except ConvergenceWarning:
            runs_res = runs_model.fit(reml=False, method="cg", start_params=runs_res.params)
        aic_dict['all_runs'].append(runs_res.aic)

    aic_df = pd.DataFrame(data=aic_dict, columns=aic_dict.keys())
    aic_df.to_csv(path_or_buf=os.path.join(TYPE_RES_DIR, type_contrast, 'cluster_AICs.csv'), index=False)

for effect, effect_contrast in zip(EFFECTS, EFFECT_CONTRASTS_OF_INTEREST):

    # dictionary saving the AIC values for each model
    aic_dict = {
        'cluster': [],
        'pre_post': [],
        'runs_in_block': [],
        'all_runs': []
    }

    # FIRST DO MAGIC BEFORE VS MAGIC AFTER
    img_path = os.path.join(EFFECT_RES_DIR, effect_contrast, CLUSTER_IMAGE)
    # read in a dataframe containing the video condition, the prepost condition and the corresponding beta names
    # created by 'create_type_ROI_DF.py'
    df_path = os.path.join(DERIVATIVES_DIR, 'univariate-ROI', 'SpecialMoment', 'MagicEffects', 'data_frame.csv')
    dataframe = pd.read_csv(df_path)
    # reduce dataframe to only Magic conditions
    magic_df = dataframe[(dataframe.Type == 'Magic') & (dataframe.Effect == effect)]

    # get all subjects
    subjects_path = os.path.join(SPM_FLA_DIR, 'MagicEffects', '6mm-smoothed-mnispace', time_analyzed)
    subjects = glob.glob(os.path.join(subjects_path, 'sub-*'))
    ALL_SUBJECTS = np.linspace(1, len(subjects), len(subjects), dtype=int)

    # create variables that are the same for all clusters
    subs = np.repeat(ALL_SUBJECTS, len(ALL_RUNS))
    runs = np.tile(ALL_RUNS, len(subjects))
    pre_post_inblock = (runs - 1) // 2  # Resulting in 0,0,1,1,2,2,3,3 ...
    # pre_post_inblock%2 resulting in 0,0,1,1,0,0,1,1 ...
    pre_post = ['pre' if pp % 2 == 0 else 'post' for pp in pre_post_inblock]
    runs_in_block = (runs - 1) % 4  # Resulting in 0,1,2,3,0,1,2,3, ...

    current_img = load_img(img_path)
    img_data = current_img.get_fdata()

    # get the highest value in our clusters image - to get the number of separate clusters
    num_clusters = np.max(img_data, axis=None).astype(int)

    for c in range(num_clusters):
        # get cluster mask
        cluster_mask = np.array(img_data == c + 1)
        # check wether the cluster is large enough
        if np.sum(cluster_mask, axis=None) < CLUSTERSIZE_THRESHOLD:
            continue

        # fill dictionary with entries
        data_dict = {
            'subject': subs,
            'runs': runs,
            'pre_post': pre_post,
            'runs_in_block': runs_in_block,
            'betas': []
        }

        for s, sub in enumerate(subjects):
            sub_df = magic_df[magic_df.subject == s + 1]
            for current_run in ALL_RUNS:
                current_beta = sub_df.BetaNames[sub_df.Runs == current_run].values[0]
                beta_nii = nib.load(os.path.join(sub, current_beta))
                beta_data = beta_nii.get_fdata()
                beta_roi = beta_data[cluster_mask & ~np.isnan(beta_data)]

                data_dict['betas'].append(beta_roi.mean(axis=None))

        cluster_df = pd.DataFrame(data=data_dict, columns=data_dict.keys())
        cluster_df.to_csv(
            path_or_buf=os.path.join(EFFECT_RES_DIR, effect_contrast, 'cluster{:02}_df.csv'.format(c + 1)), index=False)

        aic_dict['cluster'].append(c + 1)

        pre_post_model = smf.mixedlm("betas ~ pre_post", cluster_df, groups=cluster_df['subject'])
        try:
            pre_post_res = pre_post_model.fit(reml=False, method="cg")
        except ConvergenceWarning:
            pre_post_res = pre_post_model.fit(reml=False, method="cg", start_params=pre_post_res.params)
        aic_dict['pre_post'].append(pre_post_res.aic)

        runs_in_block_model = smf.mixedlm("betas ~ runs_in_block", cluster_df, groups=cluster_df['subject'])
        try:
            runs_in_block_res = runs_in_block_model.fit(reml=False, method="cg")
        except ConvergenceWarning:
            runs_in_block_res = runs_in_block_model.fit(reml=False, method="cg", start_params=runs_in_block_res.params)
        aic_dict['runs_in_block'].append(runs_in_block_res.aic)

        runs_model = smf.mixedlm("betas ~ runs", cluster_df, groups=cluster_df['subject'])
        try:
            runs_res = runs_model.fit(reml=False, method="cg")
        except ConvergenceWarning:
            runs_res = runs_model.fit(reml=False, method="cg", start_params=runs_res.params)
        aic_dict['all_runs'].append(runs_res.aic)

    aic_df = pd.DataFrame(data=aic_dict, columns=aic_dict.keys())
    aic_df.to_csv(path_or_buf=os.path.join(EFFECT_RES_DIR, effect_contrast, 'cluster_AICs.csv'), index=False)