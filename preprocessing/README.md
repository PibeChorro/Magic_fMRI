# BIDS_preprocessing
A matlab prepocessing pipeline to create a BIDS compatible data set <br>
_________________________________________________________________________________________________________________________________________________________
**IMPORTANT!!!** All the scripts here were not yet tested on a REAL BIDS-dataset and for now only works on a very simple dataset.<br>
It still has to be tested for datasets containing additional `tsv` and `json` files.<br>
_________________________________________________________________________________________________________________________________________________________

0: You need to have a project folder containing a "sourcedata" folder. This folder should contain your subject folders which themselves only contain the DICOMs<br>
- myProject
  - sourcedata
    - sub01
      - SUB001.....IMA
      - ...
      - ...
    - sub02
    - sub03

1: Run `sortDicomsIntoFolders.m` and select your "source_data" folder. The script does the following:<br>
1a: For every subject it sorts all DICOMs into folders named after the sequencenumber in which they were generated (i.e. 01 02 ...)<br>
1b: It then enters every created folder, reads out the sequence information from the first DICOM header and if this sequence is unknown (for the very first subject it will be unknown), it asks you to give it a proper name (functional data = 'func', anatomical = anat, diffusion wheight imaging = dwi, the rest you can call as you like).
You also need to deliver the number of images in this scan. ***!!!IMPORTANT!!!*** Make sure you deliver the right number of images. <br>
1c: Every folder is then moved to its corresponding sequence folder and named run-XX (if a folder contains more or less images than previously specified it will be saved as error-run-XX)<br>
1c-extra: If a sequence is called 'anat' or 'struct' you are asked to assign the modality of your anatomical scan (T1, T2 or T2star). The folder is then **RENAMED** to 'anat/<modality>' **NOT MOVED** into a 'anat/struct' folder and named 'run-XX'<br>
1d: It saves the delivered sequence information, names and coresponding number of images

- myProject
  - sourcedata
    - sequenceInfo.mat
    - sub01
      - anat
      - func
        - run-01
        - error-run-02
        - run-02
        - ...
      - loc
      - shim
    - sub02
    - sub03

2: Run `DICOMconversionSPM12_BIDS.m` (you may previously want to specify which data type you want to convert - for now only functional and anatomical data) The script does the following:<br>
2a: It askes you to select your project folder (the one CONTAINING your sourcedata folder) and to specify your subject prefix for the DICOMS.<br>
2b: The script creates a 'rawdata' folder in your project folder. This one is the folder you really want to be BIDS conform.<br>
2c: Your DICOMs in the 'sourcedata' folder will be converted into NIfTIs and stored in the 'rawdata' folder (only the 'non-error' data folders are included).<br>
2d: The created NIfTIs are then compressed into 4D NIfTIs and the existing 3D NIfTIs are deleted.<br>
2di: The 4D NIfTI files are named BIDS conform. To this point it only works for one session and one task.<br>
2e: For functional and anatomical scans a .json sidecar file is created. For now only one .json file is created per subject (we assume the same parameters for all runs).<br>

- myProject
  - rawdata
    - sub-01
      - func
        - sub-01_task-\<label\>.json
        - sub-01_task-\<label\>_run-\<index\>_bold.nii
        - sub-01_task-\<label\>_run-\<index\>_bold.nii
        - ...
      - anat
        - sub-01_\<modality\>.json
        - sub-01_\<modality\>.nii
    - sub-02
    - sub-03
  - sourcedata

**Before you can run `preprocessingSPM12_BIDS`you need to adjust the location of SPMs TMP.nii file**<br>
3: Run `preprocessingSPM12_BIDS.m` (you may previously want to specify preprocessing step you want to perform and which smoothing size you want to use) The script does the following:<br>
3a: It askes you to select your project folder (the one CONTAINING your sourcedata and rawdata folder).<br>
3b: It creates a 'derivatives' folder.<br>
3c: For every subject it creates an 'MRI' folder and within a 'func' and an 'anat' folder.
3d: For every preprocessing step it creates a corresponding folder and stores the resulting NIfTIs in their.

- myProject
  - derivative_data
    - \<other-pipeline\>
    - spm12
      - spm12-preproc
        - 6smoothed
          - sub-01
            - func
              - s6wausub-01_task\<label\>_run-\<index\>_bold.nii
              - ...
          - sub-02 
        - coregistered
        - normalized
        - realigned
        - segmented
        - slice_time_corrected
  - rawdata
  - sourcedata
  

# ROI generation
So far we have three different ROI sets, that are all generated with different atlases, toolboxes or prewritten scripts
1. Benson ROIs. These ROIs are associated with visual processing. Benson wrote a toolbox to derive these ROIs from the Wang atlas 2014
2. Glasser ROIs. These ROIs derive from the whole brain parcelation done by Glasser 2016. We chose this atlas because the parcelation is fine graned and we could use it to derive ROIs based on previous results by Danek 2015
3. Network ROIs. These networks are taken from the Yeo atlas provided by freesurfer. It gives 7 different networks
  
## 1 Benson ROIs
In order to create these ROIs, one must previously run the `recon-all` command from freesurfer. Then the `neuropythy atlas`command can be run with python. 
This script creates a two atlas volume files for each hemisphere. With `mri_cor2label` individual label files for 25 ROIs (per hemisphere) can be created. One can merge these labels if necessary/wished. Then with `mri_label2vol` one can create nifti ROI files
  
## 2 Glasser ROIs
For these ROIs its necessary to have the fsaverage image (or a link to it) in your subject directory, as well as the FreeSurferColorLUX.txt file.
Download the `create_subj_volume_parcelation.sh` file and run it for each subject. 
Depending on the flags set it creates labels and volume ROIs for each brain region in its parcelation
One then probably wants to merge some of the labels/volumes because the provided ROIs are very fine grained
  
## Network ROIs
For the seven network ROIs one needs to have access to the freesurfer directory. Freesurfer contains an atlas file (in multiple resolutions and stuff).
With the `NetworkROIs.m` script this atlas is resliced into subject nativespace and coregistered to the meanEPI image. Both images are saved in an 'atlas' folder inside the freesurfer subject directory
Running `BrainNetworkROIs.py` the atlas image is split into 7 network images, which are saved in 'corrected_ROIs'. 
*Since the original atlas was resliced, the values for the networks are not preserved, but vary. Therefore a THRESHOLD is applied assigning voxels to a network if their values are within the THRESHOLD boundaries.
  
# Event TSV files
This is not necessary for all analyses, but provides the BIDS coherent event files, containing stimulus and response information.
Run `createEventTSVfiles.m` to create the tsv files. `calculateHits.m` is also needed in that step.
