#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# Befor you run this script you have to have previously run the script 
# 'univariate_ROI-mean_DF.py'.
# The data for this analysis is derived from an fMRI experiment in which
# subjects viewed different videos. Either magic, control or surprise videos. 
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
# The purpose of this script is to read in the dataframe previously created by
# 'univariate_ROI-mean_DF.py' and based on that dataframe create a new 
# dataframe, that summarizes the data and performs some univariate analyses, 
# namely repeated measure ANOVA, on it.

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# FIRST STEP
# After getting all important path information and names of ROIs, magic EFFECTS
# (Appear, Change, Vanish), OBJECTS (Ball, Card, Stick) and video TYPES (Magic
# Control, Surprise) it creates a dictionary with all ROIs, subject, pre_post 
# and Effect as keys, and every key gets an empty list.
# SECOND STEP
# It then reads in the data frame created by 'univariate_ROI-mean_DF.py'.
# Then it iterates over subjects (outer loop), magical effects (middle loop)
# and then over ROIs (inner loop). It then takes all the data from the current
# subject, effect and ROI, that belongs to the PRE revelation condition for 
# the magic videos and calculates the mean. It does the same for the control
# videos and subtracts mean_magic - mean_control
# The result is added to the list corresponding to the current ROI in the 
# dictionary. After the inner loop the subject ID, the effect and the string
# 'pre' are added the dictionary lists.
# The same is done one more time but for the data belonging to the POST 
# revelation condition and again the subject ID, the effect and the string
# 'post' are added the dictionary lists.
# The filled dictionary is converted to a new dataframe, which is saved as a
# hdf5/csv file
# THIRD STEP
# The script iterates over all ROIs and performes a repeated measures ANOVA, 
# using the ROI data as DV and the effects and pre post condition as IV (it
# checks for sphericity). The results are saved as csv files.

######################
# COMMAND LINE FLAGS #
######################
# --analyzed, const='moment',  default='moment'
# The GLM was either calculated using all scans during one trial or depending
# on a specific moment during each video. This flag decides which of the GLMs
# is used

#############
# LIBRARIES #
#############

# interact with the operating system 
import os
import argparse
from pathlib import Path
import git
# data structuration and calculations
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
import pingouin as pg
# optimize time performance
import time

# get start time
T_START = time.time()

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--smooth",     nargs='?', const=0,         
                    default=0,          type=int)    # what data should be used
parser.add_argument("--analyzed",   nargs='?', const='moment',  
                    default='moment',   type=str)    # what data should be used
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
ANALYZED        = ARGS.analyzed

if ANALYZED == 'moment':
    data_analyzed = 'SpecialMoment'
elif ANALYZED == 'video':
    data_analyzed = 'WholeVideo'

# variables for path selection and data access
HOME                = str(Path.home())
PROJ_DIR            = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
RAWDATA_DIR         = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR     = os.path.join(PROJ_DIR, 'derivatives')
RESULTS_DIR         = os.path.join(DERIVATIVES_DIR, 'univariate-ROI',
                               data_analyzed,'VideoType')

# define ROIs
ROIS = [
        'V1', 'V2', 'V3', 'hV4', 
        'V3A', 'V3B', 
        'LO', 'VO', 
        'FEF', 'IPS',
        'ACC', 'PCC', 
        'IFG', 'aINSULA', 
        'IFJ', 'PHT', 'PF',
        'DMN', 'DAN', 'VAN', 'visual'
      ]

OBJECTS = [
    'Ball',
    'Card',
    'Stick'
    ]

VIDEO_TYPE = [
    'Magic',
    'Control',
    'Surprise'
    ]

averaged_dict = {}
for roi in ROIS:
    averaged_dict[roi] = []
    
averaged_dict['subject']     = []
averaged_dict['pre_post']   = []

# read in the previously created data frame
DATA_DF = pd.read_csv(os.path.join(RESULTS_DIR,'data_frame.csv'))
# from read data frame iterate over subjects 
# within each subject average the magic effects and the corresponding controls
# substract magic - control
# do this once for the pre revelation runs and once for the post revelation runs
# save a value for each effect pre and post in a dictionary and transform it
# into a data frame
for s,sub in enumerate(np.unique(DATA_DF.subject)):
    # get a subset of the data frame for the current subject
    sub_data_frame = DATA_DF[DATA_DF.subject==sub]
    # iterate over rois PRE revelation
    for roi in ROIS:
        mag_eff = sub_data_frame[roi][(sub_data_frame.Type=='Magic')&   # only magic data
                                      (sub_data_frame.PrePost=='pre')]      # only pre revelation
        # average the pre revelation magic effect data
        avg_mag = np.mean(mag_eff)
        
        con_eff = sub_data_frame[roi][(sub_data_frame.Type=='Control')& # only control data
                                      (sub_data_frame.PrePost=='pre')]      # only pre revelation
        # average the pre revelation controll effect data
        avg_con = np.mean(con_eff)
        
        # get difference and add it to the array of the current roi in dictionary
        mag_minus_con = avg_mag-avg_con
        averaged_dict[roi].append(mag_minus_con)
        
    # once every roi is done, add the subject ID and pre to the dict
    averaged_dict['subject'].append(s+1)
    averaged_dict['pre_post'].append('pre')
    
    # exactly the same as in the loop over ROIs above, but for post revelation data
    for roi in ROIS:
        mag_eff = sub_data_frame[roi][(sub_data_frame.Type=='Magic')&   # only magic data
                                      (sub_data_frame.PrePost=='post')]      # only post revelation
        avg_mag = np.mean(mag_eff)
        
        con_eff = sub_data_frame[roi][(sub_data_frame.Type=='Control')& # only control data
                                      (sub_data_frame.PrePost=='post')]      # only post revelation
        avg_con = np.mean(con_eff)
        
        mag_minus_con = avg_mag-avg_con
        averaged_dict[roi].append(mag_minus_con)
    averaged_dict['subject'].append(s+1)
    averaged_dict['pre_post'].append('post')
    
# convert dictionary into a pandas DataFrame for further analysis
average_df = pd.DataFrame(averaged_dict, columns=averaged_dict.keys())
average_df.to_csv(path_or_buf=os.path.join(RESULTS_DIR,'prepared_dataframe.csv'))

# since we may have to use parametric or non-parametric tests, we set up a 
# dictionary that stores the test, statistic, p-value and effect size
test_results = {
    'ROI':[],
    'test':[],
    'statistic':[],
    'p-val':[],
    'effect-size':[]
    }
# iterate over ROIS
for r,roi in enumerate(ROIS):
    # now iterate over the rois you want to use for the analysis
    test_results['ROI'].append(roi)
    x_var = average_df[average_df.pre_post=='pre'][roi]
    y_var = average_df[average_df.pre_post=='post'][roi]
    
    norm = pg.normality(data=average_df,
                        dv=roi,
                        group='pre_post')
    if all(norm.normal.values):
        # parametric t-test
        ttest_res = pg.ttest(x=x_var,y=y_var,paired=True)
        test_results['test'].append('t-test')
        test_results['statistic'].append(ttest_res['T'].values[0])
        test_results['p-val'].append(ttest_res['p-val'].values[0])
        test_results['effect-size'].append(ttest_res['cohen-d'].values[0])
    else:
        # do non-parametric test
        wilcox_res = pg.wilcoxon(x=x_var, y=y_var)
        test_results['test'].append('wilcoxon')
        test_results['statistic'].append(wilcox_res['W-val'].values[0])
        test_results['p-val'].append(wilcox_res['p-val'].values[0])
        test_results['effect-size'].append(wilcox_res['RBC'].values[0])
    
results_df = pd.DataFrame(data=test_results,columns=test_results.keys())
results_df['p-corrected'] = results_df['p-val']*len(results_df)
results_df['p-corrected'][results_df['p-corrected']>1] = 1

results_df.to_csv(path_or_buf=os.path.join(RESULTS_DIR,'test_results.csv'),
                  index=False, float_format='%.3f')
        
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
    
# create a log file, that saves some information about the run script
with open(os.path.join(RESULTS_DIR,'analysis-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))
