function smoothFmriprep()

%% Preprocess fMRI data using spm12
% JB 08/2019 (adapted from PRG 05/2019)
% Further changed by VP 01/2021

%% Define important details of your file structure and location
homedir = '/';
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "raw_data")\n\n'])
rootDir    = uigetdir(homedir, 'Select Project Folder');
if rootDir == 0
    error('No folder was selected --> I terminate the script')
end

% Set source_data directory. This is only needed for slice time correction
sourceDir  = fullfile(rootDir,'sourcedata');
if ~isfolder(sourceDir)
    fprintf(['It appears you do not have a "source_data" folder.\n'...
        'Please select the folder that contains your unprocessed niftis.\n\n'])
    sourceDir  = uigetdir(rootDir, 'Select DICOM folder');
    if sourceDir == 0
        error('No folder was selected --> I terminate the script')
    end
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

derivativesDir = fullfile (rootDir, 'derivatives');
if ~isfolder(derivativesDir)
    mkdir (derivativesDir);
end

subPrefix      = 'sub-';   % subject prefix - this one is REQUIRED
sessionPrefix  = 'ses';    % session prefix - only needed if more than one session per subject
sessionLabels  = '';       % labels of sessions - must be alphanumeric
taskPrefix     = 'task-';  % this one is REQUIRED (it does not make too much sense, when having only one task, but it gives more information in the name
taskLabels     = {'magic'};  % labels of tasks
anatDir        = 'anat';   % where the structural image is stored
anatModality   = 'T1w';

%% Set name of data
% A first draft of how one could automatise the names of data
% sessionNames = {};
% if isempty(sessionLabels)
%     sessionNames = '';
% else
%     for ses = 1:length(sessionLabels)
%         sessionNames{end+1} = [sessionPrefix sessionLabels{ses}];
%     end
% end
% 
% taskNames = {};
% for task = 1:length(taskLabels)
%     taskNames{end+1} = [taskPrefix taskLabels{task}];
% end
% 
% dataNames = {};
% for ses = 1:length(sessionNames)
%     for task = 1:length(taskNames)
%         dataNames{end+1} = [sessionNames{sees} '_' taskNames{task}];
%     end
% end

taskName = [taskPrefix taskLabels{1}];

%-------------------------------------------------------------------------%


%% Decide what to do
%..............................WHAT TO DO.................................%
do.overwrite        = 1;
do.gunzip           = 1;
do.smoothing        = 1; %Smoothing Flag, set to 1 if you want to smooth. 
do.smoothNorm       = 'mni'; % Smooth normalize data = 'mni', native data = 'native' or 'both'
do.smoothingSize    = 6; % in mm 

%% OPEN SPM -- necessary ?
spm fmri;

% create a BIDS conform directory structure for the NIFTIS
% first we need to create a cell containing the subject names
softwareName    = 'fmriprep-21.0.1';
preprocessedDir = fullfile(derivativesDir,softwareName);
subNames        = cellstr(spm_select('List', preprocessedDir, 'dir',['^' subPrefix]));

%% start to perform the preprocessing
for ss = 2:length(subNames) % For all subjects do each ...


    currDir = fullfile(preprocessedDir, subNames{ss}, 'func');

    if do.gunzip
        runs = cellstr(spm_select('FPList', currDir, 'space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')); 
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.files = runs;
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.outdir = {currDir};
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.keep = true;
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear matlabbatch
    end
    
    %% Smoothing
    if do.smoothing
        fprintf('SMOOTHING\n\n')
        if strcmp(do.smoothNorm,'mni') == 1
            runs                = cellstr(spm_select('List', currDir, 'space-MNI152NLin2009cAsym_desc-preproc_bold.nii$'));
            runs                = natsort(runs);
        elseif strcmp(do.smoothNorm, 'native') == 1
            folderContent = dir(fullfile(currDir, subNames{ss}, 'func', ['au' subNames{ss} '_' taskName '_' 'run*.nii']));
        end

        nruns         = length(runs);

        alltargets = {};

        for r = 1:nruns
            dirfiles     = spm_select('ExtFPList', currDir, runs{r}, Inf);

            if strcmp(dirfiles,'')
                warning('No files selected!');
                return;
            end
            alltargets = [alltargets; cellstr(dirfiles)];
        end

        matlabbatch{1}.spm.spatial.smooth.data      = cellstr(alltargets);
        matlabbatch{1}.spm.spatial.smooth.fwhm      = repmat(do.smoothingSize,1,3); % Set at the beginning.
        matlabbatch{1}.spm.spatial.smooth.dtype     = 0;
        matlabbatch{1}.spm.spatial.smooth.im        = 1;
        matlabbatch{1}.spm.spatial.smooth.prefix    = ['s' num2str(do.smoothingSize)];

        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clearvars matlabbatch
%         [success,message] = movefile(fullfile(currDir, ['s' num2str(do.smoothingSize) '*.nii']),...
%             fullfile(smoothedDir, subNames{ss}, 'func'));
%         if ~success
%             warning(message)
%         end
    end
end
end