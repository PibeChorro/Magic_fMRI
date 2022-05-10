# Magic_fMRI
all around the Magic fMRI project started 2020

# General
Here EVERYTHING that has to do ANYTHING with the Magic fMRI project startet 2020 SHOULD be included.

There is NO clear linearity in all the scripts and some analyses are independent from one another.<br>
However, some steps need to be done in the right order. Normally scripts in the ```preprocessing``` folder need to be run first. Only then, scrips from the ```analysis``` folder can be run

## preprocessing
### Very first steps
In the beginning MRI DICOM files **must** be located in ```sourcedata``` and subjects foulders are called ```sMag<index>``` (better would have been a BIDS naming). 

DICOM files are of the form ```SMAG<index>_<date>.MR.MPI_STANDARD_PABPLO.<sequence-nr>.<volume-nr>.<year>.<month>.<day>.<series of numbers(meaning unknown)>.IMA```

```sortDicomsIntoFolders.m``` asks for the directory of your sourcedata (where DICOMS are located), the prefix of your subjects (sMag in this case) and for the number of dummy scans, that should be excluded. Based on the sequence number the data is moved into seperate folders. <br>
After sorting DICOMS into sequence folders, the script iterates over the folders, reads out the *name* of the sequence from the first DICOM header and asks you to give data from this sequence an apropriate name and how many volumes these sequence should have measured. It then renames the folders accordingly

After sorting DICOMS into folders ```DICOMconversionSPM12_BIDS``` uses SPM12 to convert DICOMS into NIfTIs, move them into a newly created ```rawdata``` folder and compresses the multiple 3D files into a single 4D file.

**MAKE SURE THAT SORTING AND CONVERSION WERE PERFORMED CORRECTLY**

### Further steps
The next few steps are somewhat independent from each other.

1) SPM preprocessing

For wholebrain analysis of the data ```preprocessingSPM12_BIDS.m``` should be run. It runs the typical preprocessing steps in the following order
- realign the NIfTIs from rawdata
- perform slicetime correction from realigned data
- coregister slicetime corrected data to the anatomical image
- segmentation -> split the anatomical image into a graymatter, whitematter and CSF mask
- normalize coregistered data into MNI space
- smooth data with specified smoothing kernel (6mm in this case)

For each of the steps the created data is moved into a seperate folder to prevent further analyses from using the wrong data.

2) fmriprep preprocessing (**not in use right now**)<br>
```preprocessing_fmriprep.sh```
This file runs a docker image (that can also be run on the cluster) to preprocess the fMRI data. fmirprep performes more steps, but does not create seperate files for each step. Different flags can be set in order to add or surpress some steps (IMPORTANT!!! right now the script DOES NOT perform the freesurfer reconstruction)

3) freesurfer reconstruction and ROI generation

```docker_fs.sh```<br>
This file also runs a docker freesurfer image. It runs some different steps:
- create a directory for each subject and copy the anatomical image into it
- perform the ```recon-all -autorecon-all``` comand (takes around a day)
- run ```neuropythy atlas``` (not a freesurfer, but a python command - https://github.com/noahbenson/neuropythy)
- the created mgz file is transformed into label files (for each hemisphere)
- (ditched for the moment) create a ```registration.dat``` file with ```bbregister``` for registering a surface label file to a volume  
- label files are transformed into volume files (**IMPORTANT! the transformation need an image for a reference space - the mean coregistered image created by SPM is used**)
- additionally create a volume for the 3rd ventricle as a control ROI

```GlasserROIs.sh```
Based on the results from Danek et al and Parris et al we wanted to add ROIs to our set.<br>
This file creates a docker freesurfer container that merges selected labels from the HCP Glasser atlas and transforms these into volumes

3a) *Additional step:* After creating ROIs, run ```assign_ROI_overlaps.py``` to make sure each voxel only belongs to one ROI.

4) network ROIs

How the different networks in the brain (DAN, VAN, DMN and visual network) react to different conditions was of interest
```NetworkROIs.m```
This file takes the Yeo atlas from freesurfer directory and transforms it into functional space of each subject
```BrainNetworkROIs.py```
Since the volumetric Yeo atlas sets a value to each voxel to assign it to a network and this volume is transformed, the values are also changed from easily distinctable integers to floats. This script reassignes the values based on thresholds 

5) creation of event TSV files

```createEventTSVfiles.m```
This file goes through the log files created by the experiment and creates event tsv file for each run.<br>
For specific subjects missing values are replaced with ```4``` since the button did not work properly in some measurements.

## analysis
The majority of scripts here are independent from each other. However, some scripts create files that are needed in other scripts (e.g., data frame like csv files)<br>

### Whole Brain analysis
When ```preprocessingSPM12_BIDS.m``` has been run one can run 
- ```spm_fla_VideoTypes_BIDS.m``` estimates each video type (magic, control and surprise)
- ```spm_fla_MagicEffects.m``` estimates each magic effect for the magic and control videos (appear, change, vanish) and surprise videos
- ```spm_fla_EveryVideo.m``` estimates each video on its own
Each runs a first level analysis for each subjecect, with response times, movement parameters and WM/CSF/GS as regressors of no interest<br>
IMPORTANT! there are different branches with different regressors of no interest. Which one is supposed to be the right is not yet decided<br>
The first two create contrast images that can later be used for second level analysis

When ```spm_fla_VideoTypes_BIDS.m``` and /or ```spm_fla_MagicEffects.m``` has been run
- ```spm12sla_BIDS.m``` performs aparametric second level analysis on all contrasts calculated in the first level analysis
- ```snpm13sla.m``` performs a non-parametric second level analysis by means of a sign-flipping test on all contrasts calculated in the first level analysis

When ```spm12sla_BIDS.m``` and /or ```snpm13sla.m``` has been run on the Effect specific first level contrasts
- ```Conjunction_Effects.m```  reads in T-value maps (parametric) or -log10 p-value maps (nonparametric) and looks for overlaps in the effects -> decision for parametric or non-parametric and corresponding values has to be made

When ```spm_fla_EveryVideo.m``` has been run in MNI space
- ```decode_magic_effects_searchlight.py``` runs a searchlight decoding analysis in a single subject using LDA to decode the magic effect
- ```decode_magic_vs_nomagic_searchlight.py``` runs a searchlight decoding analysis in a single subject using LDA to decode if magic or no-magic was shown<br>
both scripts might need some fine tuning.

When ```decode_magic_effects_searchlight.py``` and /or ```decode_magic_vs_nomagic_searchlight.py``` has been run<br>
IMPORTANT! this is an order!
- ```searchlight_bootstrap_procedure.py``` creates bootstrapped NIfTI images (10000)
- ```searchlight_cluster_statistic.py``` uses the bootstrapped images to calculate a non-parametric corrected p-value for searchlight decoding accuracies (**not correct! Needs to be checked again**)<br>
Also ```searchlight_parametric_group_analysis.py``` can be run, but is less accurate

### ROI analysis
**ALL** scripts in this subsection require at least one of the SPM first level analysis results (in native space)

#### univariate analysis
There are three python scripts that generate dataframes with values for each ROI based on first level beta estimates
- ```create_type_ROI_DF.py```
- ```create_effect_ROI_DF.py```
- ```create_Rating_Betas_DF.py```
For each of the estimated beta (of interest) the script calculates the mean value of voxels within a ROI and add some information (e.g., run, pre post condition etc.)<br>
```create_Rating_Betas_DF.py``` also adds the subjective surprise rating for each trial

When ```create_type_ROI_DF.py``` and ```create_effect_ROI_DF.py``` has been run, the corresponding ```type_ROI_analysis.py``` and ```effect_ROI_analysis.py``` can be run.<br>
When ```create_Rating_Betas_DF.py``` has been run ```rating_ROI_correlation.py``` can be run (the exact correlation has to be discussed again!)<br>
ALSO ALSO after ```create_Rating_Betas_DF.py``` the behavioral analysis can be run with ```Rating_analysis.py``` (maybe I'll change ANOVAS and Friedman Tests to mixed effects models)

When ```spm12sla_BIDS.m``` and /or ```snpm13sla.m``` has been run 
- ```create_ROIs_from_results.py``` creates ROIs from significant clusters in the prepost and interaction condition (located in ```preprocessing```)
- ```timeevolution_clusters.py``` checks weather the significant results stem from the difference in prior knowledge or from a time confound

#### decoding analysis
When ```spm_fla_EveryVideo.m``` has been run
- ```decode_magic_effects.py```
- ```decode_magic_pre_vs_post.py``` (not yet decided how exactly)
- ```decode_magic_vs_nomagic.py```

When one of the above decoding scripts has been run
- ```permutation_max_statistic.py``` (split ROI groups)

## eyetracking
Before eyetracking can be performed, the ```.edf``` file created by the eyetracker has to be converted into an ```asc``` file. This can be performed with the SR-Research ```edf2asc``` script.

First step:
preprocess ET data with ```preprocessing_ET_data.py``` this does the following
- remove data before and after blinks
- interpolate missing data
- low pass and high pass filter pupil diameter data
- memean pupil diameter per condition
- divide by standard deviation per session (run)
- downsample all data to video framerate

After that ```saccades_blinks_df.py``` can be run. Then use the created dataframe to run
- ```blink_saccade_analysis.py```
- ```compare_prepost_fixations.py```
# TODO
Add edf2ascii to the code.
