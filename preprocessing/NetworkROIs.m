%% HEADER
% The purpose of this script is to transform the Yeo 7 Network atlas into
% native space for each subject.
% It takes the
% Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii image
% located in the freesurfer directory under
% 'average/Yeo_JNeurophysiol11_MNI152' (you need to unzip the image first)
% and the inverse normalization image produced in the segmentation
% preprocessing step performed by spm.
% For each subject after transformation the image is moved into the
% freesurfer directory of each subject (in an 'atlases' folder, which is
% created if not existing)

%% START
clear all
close all

%% Define important details of your file structure and location
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "raw_data")\n\n'])
rootDir    = '/Users/vpl/Documents/Master_Thesis/DATA/MRI'; % uigetdir(homedir, 'Select Project Folder');
if rootDir == 0
    error('No folder was selected --> I terminate the script')
end

% Set rawdata directory.
rawDir     = fullfile(rootDir, 'rawdata');
if ~isfolder(rawDir)
    fprintf(['It appears you do not have a "rawdata" folder.\n'...
        'Please select the folder that contains your unprocessed niftis.'])
    rawDir  = uigetdir(rootDir, 'Select unprocessed nifti folder');
    if rawDir == 0
        error('No folder was selected --> I terminate the script')
    end
end

% Set derivatives directory
derivativesDir = fullfile (rootDir, 'derivatives');
if ~isfolder(derivativesDir)
    mkdir (derivativesDir);
end

% Get the full path directories for each subject of the segmented
% preprocessing step and the names of the subjects for later data structure
segmentedDir    = fullfile(derivativesDir,'spm12', 'spm12-preproc','segmented');
segmentedSubs   = cellstr(spm_select('FPList',segmentedDir,'dir','sub-'));

coregisteredDir = fullfile(derivativesDir, 'spm12', 'spm12-preproc','coregistered');
coregisteredSubs= cellstr(spm_select('FPList',coregisteredDir,'dir','sub-'));

subNames        = cellstr(spm_select('List',segmentedDir,'dir','sub-'));

% Get the path where the atlas is located
pipelineName    = 'freesurfer';
atlasDir        = fullfile(homedir,'Documents',pipelineName, 'average', 'Yeo_JNeurophysiol11_MNI152');
yeoAtlas        = 'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii';

spm('defaults', 'FMRI');

%% Iterate over subjects and transform the atlas into subject nativespace
for s = 1:length(subNames)
    % Where to store/move the native space altas
    destDir = fullfile(derivativesDir,pipelineName,subNames{s},'atlases');
    if ~isfolder(destDir)
        mkdir (destDir);
    end
    % get the deformation image
    defFieldDir = fullfile(segmentedSubs{s}, 'anat');
    defField    = spm_select('FPList', defFieldDir,'iy_sub-');
    
    % normalization settings (all default, except for 'vox')
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(fullfile(atlasDir,yeoAtlas));
    % set deformation image in the normalization settings
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(defField);
    
    spm_jobman('run', matlabbatch);
    
    % move created image to subject directory
    [success,message] = movefile(string(fullfile(atlasDir, ['w' yeoAtlas])), ...
            fullfile(destDir, 'Yeo_7Network.nii'),'f');
    if ~success
        warning(message)
    end
    clearvars matlabbatch
    
    % coregister to mean EPI image and reslice to get same dimensions
    % get mean EPI image
    meanEPIDir  = fullfile(coregisteredSubs{s}, 'func');
    meanEPI     = spm_select('FPList', meanEPIDir, 'mean');
    
    % coregister specifications
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(meanEPI);
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(fullfile(destDir, 'Yeo_7Network.nii'));
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    
    spm_jobman('run', matlabbatch);
    clearvars matlabbatch
end

%% END
clear all
close all