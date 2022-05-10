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

2) fmriprep preprocessing (**not in use right now**)
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
Most simple are the univariate whole brain analyses and only need the ```preprocessingSPM12_BIDS``` script run
decoding and ROI analyses need freesurfer stuff done

# TODO
Add edf2ascii to the code.
