# Magic_fMRI
all around the Magic fMRI project started 2020

# General
Here EVERYTHING that has to do ANYTHING with the Magic fMRI project startet 2020 SHOULD be included.

There is NO clear linearity in all the scripts and some analyses are independent from one another.
However, some steps need to be done in the right order. Normally scripts in the ´´´preprocessing´´´ folder need to be run first. Only then, scrips from the ´´´analysis´´´ folder can be run

## preprocessing
In any case ´´´sortDicomsIntoFolders.m´´´ and ´´´DICOMconversionSPM12_BIDS´´´ need to be run to get nifti files.
Once DICOMS are converted into nifti, the "real" work starts. First thing should be most of the stuff in ´´´preprocessingSPM12_BIDS´´´ to get *at least* realigned, slicetime corrected and coregistered data.
After that, it makes sense to run all freesurfer scripts, in order to get surface volumes, labels and potentially ROIs

## analysis
Most simple are the univariate whole brain analyses and only need the ´´´preprocessingSPM12_BIDS´´´ script run
decoding and ROI analyses need freesurfer stuff done

# TODO
Add edf2ascii to the code.
