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
# 'create_effect_ROI_DF.py' and based on that dataframe create a new 
# dataframe, that summarizes the data and performs some univariate analyses, 
# on each ROI. If normality is given do repeated measures, else Friedman and
# Wilcoxon tests

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
# .csv file
# THIRD STEP
# The script iterates over all ROIs and performes a repeated measures ANOVA, 
# using the ROI data as DV and the effects and pre post condition as IV (it
# checks for sphericity) or a Friedman test for the effects condition and a 
# Wilcoxon test for the pre-post condition if normality is not given. 
# The results are saved as csv files.

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
# library for neuroimaging
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
                               data_analyzed,'MagicEffects')

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

EFFECTS = [
    'Appear',
    'Change',
    'Vanish'
    ]

averaged_dict = {}
for roi in ROIS:
    averaged_dict[roi] = []
    
averaged_dict['subject']     = []
averaged_dict['pre_post']   = []
averaged_dict['Effect']     = []

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
    # iterate over magic effects 
    for effect in EFFECTS:
        # iterate over rois PRE revelation
        for roi in ROIS:
            mag_eff = sub_data_frame[roi][(sub_data_frame.Effect==effect)&  # only data from the current effect
                                          (sub_data_frame.Type=='Magic')&   # only magic data
                                          (sub_data_frame.PrePost=='pre')]      # only pre revelation
            # average the pre revelation magic effect data
            avg_mag = np.mean(mag_eff)
            
            con_eff = sub_data_frame[roi][(sub_data_frame.Effect==effect)&  # only data from the current effect
                                          (sub_data_frame.Type=='Control')& # only control data
                                          (sub_data_frame.PrePost=='pre')]      # only pre revelation
            # average the pre revelation controll effect data
            avg_con = np.mean(con_eff)
            
            # get difference and add it to the array of the current roi in dictionary
            mag_minus_con = avg_mag-avg_con
            averaged_dict[roi].append(mag_minus_con)
            
        # once every roi is done, add the subject ID and pre to the dict
        averaged_dict['subject'].append(s+1)
        averaged_dict['pre_post'].append('pre')
        averaged_dict['Effect'].append(effect)
        
        # exactly the same as in the loop over ROIs above, but for post revelation data
        for roi in ROIS:
            mag_eff = sub_data_frame[roi][(sub_data_frame.Effect==effect)&  # only data from the current effect
                                          (sub_data_frame.Type=='Magic')&   # only magic data
                                          (sub_data_frame.PrePost=='post')]      # only post revelation
            avg_mag = np.mean(mag_eff)
            
            con_eff = sub_data_frame[roi][(sub_data_frame.Effect==effect)&  # only data from the current effect
                                          (sub_data_frame.Type=='Control')& # only control data
                                          (sub_data_frame.PrePost=='post')]      # only post revelation
            avg_con = np.mean(con_eff)
            
            mag_minus_con = avg_mag-avg_con
            averaged_dict[roi].append(mag_minus_con)
        averaged_dict['subject'].append(s+1)
        averaged_dict['pre_post'].append('post')
        averaged_dict['Effect'].append(effect)
    
# convert dictionary into a pandas DataFrame for further analysis
average_df = pd.DataFrame(averaged_dict, columns=averaged_dict.keys())
average_df.to_csv(path_or_buf=os.path.join(RESULTS_DIR,'prepared_dataframe.csv'))

# create a effect and a pre_post main effect and interaction dictionary
effect_maineffect_results = {
    'ROI':[],
    'test':[],
    'statistic':[],
    'ddof1':[],
    'ddof2':[],
    'p-val':[],
    'effect-size':[]
    }

prepost_maineffect_results = {
    'ROI':[],
    'test':[],
    'statistic':[],
    'ddof1':[],
    'ddof2':[],
    'p-val':[],
    'effect-size':[]
    }

interaction_results = {
    'ROI':[],
    'test':[],
    'statistic':[],
    'ddof1':[],
    'ddof2':[],
    'p-val':[],
    'effect-size':[]
    }
# iterate over ROIS
for r,roi in enumerate(ROIS):
    # now iterate over the rois you want to use for the analysis
    
    # append the ROI to the effect and revelation main effect dictionaries
    effect_maineffect_results['ROI'].append(roi)
    prepost_maineffect_results['ROI'].append(roi)
    
    # test for normality in the revelation condition ...
    pre_post_norm = pg.normality(data=average_df,
                                 dv=roi,
                                 group='pre_post')
    # ... and the effect condition
    effect_norm = pg.normality(data=average_df,
                               dv=roi,
                               group='Effect')
    
    # if normality is given ...
    if all(pre_post_norm.normal.values) and all (effect_norm.normal.values):
        # ... do a repeated measures ANOVA (rm_ANOVA)
        #############################################################################
        # CAUTION: at the moment the data frame returned by pg.rm_anova is the same #
        # no matter if you provide detailed or correction                           #
        #############################################################################
        # apply a sphericity test. If the sphericity is violated use correction
        # in ANOVA
        spher = pg.sphericity(data=average_df, 
                              dv=roi, 
                              within=['pre_post','Effect'],
                              subject='subject')
        
        # add the type of test to the dictionaries
        effect_maineffect_results['test'].append('rmANOVA')
        prepost_maineffect_results['test'].append('rmANOVA')
        # add ROI to interaction dictionary
        interaction_results['ROI'].append(roi)
        interaction_results['test'].append('rmANOVA')
        
        # pingouin's rm ANOVA function
        aov_res = pg.rm_anova(data=average_df, 
                              dv=roi, 
                              within=['Effect','pre_post'],
                              subject='subject',
                              correction = not spher.spher,
                              detailed=True)
        
        # add F-value to dictionaries
        effect_maineffect_results['statistic'].append(aov_res['F'].values[0])
        prepost_maineffect_results['statistic'].append(aov_res['F'].values[1])
        interaction_results['statistic'].append(aov_res['F'].values[2])
        
        # add degrees of freedom to dicionaries
        effect_maineffect_results['ddof1'].append(aov_res['ddof1'].values[0])
        prepost_maineffect_results['ddof1'].append(aov_res['ddof1'].values[1])
        interaction_results['ddof1'].append(aov_res['ddof1'].values[2])
        
        effect_maineffect_results['ddof2'].append(aov_res['ddof2'].values[0])
        prepost_maineffect_results['ddof2'].append(aov_res['ddof2'].values[1])
        interaction_results['ddof2'].append(aov_res['ddof2'].values[2])
        
        # add p-values to dicionaries
        effect_maineffect_results['p-val'].append(aov_res['p-GG-corr'].values[0])
        prepost_maineffect_results['p-val'].append(aov_res['p-GG-corr'].values[1])
        interaction_results['p-val'].append(aov_res['p-GG-corr'].values[2])
        
        # add effect size (partiel eta square - np2) to dicionary
        effect_maineffect_results['effect-size'].append(aov_res['np2'].values[0])
        prepost_maineffect_results['effect-size'].append(aov_res['np2'].values[1])
        interaction_results['effect-size'].append(aov_res['np2'].values[2])
        
        # if we have a significant interaction we need to do post-hoc t-tests
        if aov_res['p-GG-corr'].values[2]<=.05:
            post_hoc_ttest = pg.pairwise_ttests(data=average_df,
                                                dv=roi,
                                                within=['pre_post','Effect'],
                                                subject='subject',
                                                effsize='cohen',
                                                padjust='holm')
            post_hoc_ttest.to_csv(path_or_buf=os.path.join(RESULTS_DIR,roi+'_interaction_post_hoc.csv'),
                                  index=False,float_format='%.3f')
            
        # if there is no interaction but a main effect ...
        elif aov_res['p-GG-corr'].values[0]<=.05 or aov_res['p-GG-corr'].values[1]<=.05:
            # ... for the effect, do post-hoc t-test for effects and save result as .csv
            if aov_res['p-GG-corr'].values[0]<=.05:
                post_hoc_ttest = pg.pairwise_ttests(data=average_df,
                                                    dv=roi,
                                                    within='Effect',
                                                    subject='subject',
                                                    effsize='cohen',
                                                    padjust='holm')
                post_hoc_ttest.to_csv(path_or_buf=os.path.join(RESULTS_DIR,roi+'_effect_post_hoc.csv'),
                                      index=False,float_format='%.3f')
            # ... for the revelation, do post-hoc test for revelation and save result as .csv
            if aov_res['p-GG-corr'].values[1]<=.05:
                post_hoc_ttest = pg.pairwise_ttests(data=average_df,
                                                    dv=roi,
                                                    within='pre_post',
                                                    subject='subject',
                                                    effsize='cohen',
                                                    padjust='holm')
                post_hoc_ttest.to_csv(path_or_buf=os.path.join(RESULTS_DIR,roi+'_prepost_post_hoc.csv'),
                                      index=False,float_format='%.3f')
                
    # else do non-parametric tests
    else:
        # non parametric friedman tests
        
        # add type of test to dictionaries
        effect_maineffect_results['test'].append('Friedman')
        prepost_maineffect_results['test'].append('Wilcoxon')
        
        # since there is no (easy) alternative to a two way rmANOVA, we do 
        # friedman test for the effect main effect
        friedman_res = pg.friedman(data=average_df,
                                   dv=roi,
                                   within='Effect',
                                   subject='subject',
                                   method='f')
        
        # add F-value, degrees of freedom, p-value to dictionary
        effect_maineffect_results['statistic'].append(friedman_res['F'].values[0])
        effect_maineffect_results['ddof1'].append(friedman_res['ddof1'].values[0])
        effect_maineffect_results['ddof2'].append(friedman_res['ddof2'].values[0])
        effect_maineffect_results['p-val'].append(friedman_res['p-unc'].values[0])
        effect_maineffect_results['effect-size'].append(friedman_res['W'].values[0])
        
        x_var = average_df[average_df.pre_post=='pre'][roi]
        y_var = average_df[average_df.pre_post=='post'][roi]
        roi_wilcox = pg.wilcoxon(x=x_var, y=y_var)
        
        # add W-value, degrees of freedom, p-value and effect size to dictionary
        prepost_maineffect_results['statistic'].append(roi_wilcox['W-val'].values[0])
        prepost_maineffect_results['ddof1'].append('')
        prepost_maineffect_results['ddof2'].append('')
        prepost_maineffect_results['p-val'].append(roi_wilcox['p-val'].values[0])
        prepost_maineffect_results['effect-size'].append(roi_wilcox['RBC'].values[0])
        
        # Interaction - for every non-parametric test we do "post-hoc" 
        # interaction tests
        # empty lists to store what is the x and y variable in which data to be tested
        x_var   = []
        y_var   = []
        data    = []
        # empty data frame with wilcoxon columns
        non_para_interaction_df = pd.DataFrame(columns=roi_wilcox.keys())
        
        # wilcoxon test for pre post in all effects
        for eff in EFFECTS:
            x = average_df[(average_df.pre_post=='pre') & (average_df.Effect==eff)][roi]
            y = average_df[(average_df.pre_post=='post') & (average_df.Effect==eff)][roi]
            data.append(eff)
            x_var.append('pre')
            y_var.append('post')
            wilcox = pg.wilcoxon(x=x,y=y)
            non_para_interaction_df = pd.concat([non_para_interaction_df,wilcox])
            
        # wilcoxon test for all effects in pre post 
        for pp in ['pre', 'post']:
            app = average_df[(average_df.pre_post==pp) & (average_df.Effect=='Appear')][roi]
            cha = average_df[(average_df.pre_post==pp) & (average_df.Effect=='Change')][roi]
            van = average_df[(average_df.pre_post==pp) & (average_df.Effect=='Vanish')][roi]
            
            wilcox = pg.wilcoxon(x=app,y=cha)
            data.append(pp)
            x_var.append('Appear')
            y_var.append('Change')
            non_para_interaction_df = pd.concat([non_para_interaction_df,wilcox])
            
            wilcox = pg.wilcoxon(x=app,y=van)
            data.append(pp)
            x_var.append('Appear')
            y_var.append('Vanish')
            non_para_interaction_df = pd.concat([non_para_interaction_df,wilcox])
            
            wilcox = pg.wilcoxon(x=cha,y=van)
            data.append(pp)
            x_var.append('Change')
            y_var.append('Vanish')
            non_para_interaction_df = pd.concat([non_para_interaction_df,wilcox])
            
        # add the x and y-variable as well as data to the data frame and perform
        # holmeroni correction for multiple comparisons. Save results as .csv
        non_para_interaction_df['X']=x_var
        non_para_interaction_df['Y']=y_var
        non_para_interaction_df['data']=data
        non_para_interaction_df['p-corrected']=non_para_interaction_df['p-val']*len(non_para_interaction_df)
        non_para_interaction_df['p-corrected'][non_para_interaction_df['p-corrected']>1] = 1
        
        non_para_interaction_df.to_csv(path_or_buf=os.path.join(RESULTS_DIR,roi +'_wilcox_interaction.csv'),
                      index=False, float_format='%.3f')
        
        
# turn dictionaries into dataframes, add holmeroni correction and save as .csv
pre_post_maineffect_df= pd.DataFrame(data=prepost_maineffect_results,
                                     columns=prepost_maineffect_results.keys())
pre_post_maineffect_df['p-corrected'] = pre_post_maineffect_df['p-val']*len(prepost_maineffect_results)
pre_post_maineffect_df['p-corrected'][pre_post_maineffect_df['p-corrected']>1] = 1

pre_post_maineffect_df.to_csv(path_or_buf=os.path.join(RESULTS_DIR,'pre_post_maineffect.csv'),
                              index=False, float_format='%.3f')

effect_maineffect_df= pd.DataFrame(data=effect_maineffect_results,
                                     columns=effect_maineffect_results.keys())
effect_maineffect_df['p-corrected'] = effect_maineffect_df['p-val']*len(prepost_maineffect_results)
effect_maineffect_df['p-corrected'][effect_maineffect_df['p-corrected']>1] = 1

effect_maineffect_df.to_csv(path_or_buf=os.path.join(RESULTS_DIR,'effect_maineffect.csv'),
                            index=False, float_format='%.3f')

interaction_df= pd.DataFrame(data=interaction_results,
                                     columns=interaction_results.keys())
interaction_df['p-corrected'] = interaction_df['p-val']*len(interaction_df)
interaction_df['p-corrected'][interaction_df['p-corrected']>1] = 1

interaction_df.to_csv(path_or_buf=os.path.join(RESULTS_DIR,'interaction.csv'),
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
with open(os.path.join(RESULTS_DIR,'rmANOVAs-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))
