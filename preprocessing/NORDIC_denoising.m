%% Define important details of your file structure and location
homedir = '/home/vplikat';
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "rawdata")\n\n'])
rootDir    =  uigetdir(homedir, 'Select Project Folder');
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

derivesDir = fullfile (rootDir, 'derivatives');
if ~isfolder(derivesDir)
    fprintf(['It appears you do not have a "derivatives" folder.\n'...
        'Please select the folder that contains your preprocessed niftis.'])
    derivesDir  = uigetdir(rootDir, 'Select preprocessed nifti folder');
    if derivesDir == 0
        error('No folder was selected --> I terminate the script')
    end
end

%% Data locations 
softwareName        = 'nordic';              % software used to create preprocessed data
destDir             = fullfile(derivesDir,softwareName);

subNames            = spm_select('List', rawDir, 'dir', 'sub-');
subNames            = cellstr(subNames);

ARG.magnitude_only = 1;
for s = 1
    rawSubDir       = fullfile(rawDir,subNames{s});
    rawSubFuncDir   = fullfile(rawSubDir,'func');
    outputDir       = fullfile(destDir,subNames{s});
    runs            = cellstr(spm_select('List', rawSubFuncDir, '_bold.nii')); 
    numRuns         = length(runs); 
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    
    cd(outputDir)

    for r = 1:numRuns
        input_file = fullfile(rawSubFuncDir, runs{r});
        output_file = fullfile(outputDir, runs{r});
        
        NIFTI_NORDIC(input_file,'', runs{r}(1:end-4), ARG);
    end
end