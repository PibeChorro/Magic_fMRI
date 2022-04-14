%%%%%%%%%%
% HEADER %
%%%%%%%%%%
%% Purpose of this scirpt
% Create a conjunction map, showing only significant voxels in all three
% effects (Appear, Change, Vanish) pre revelation.

%% Functionality of this script
% FIRST STEP
% Set up the necessary directories for input and output
% Select the images to be used in the conjunction analysis
% Run conjunction using SPM12

clear all; % first clear all variables, so we don't have any intervening variables.
tic; % start script.

%% Define important details of your file structure and location
homedir = '/home/vplikat';
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "rawdata")\n\n'])
rootDir    = uigetdir(homedir, 'Select Project Folder');
if rootDir == 0
    error('No folder was selected --> I terminate the script')
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
softwareName        = 'snpm13';              % software used to create preprocessed data

% specify the name of the analysis pipeline
analysisPipeline    = 'snpm13-sla';          % how is the folder named that contains first level results
brainMask           = 'WholeBrain';         % whole brain or ROI
glmUsed             = 'MagicEffects';
smoothKernelSize    = 6;                    % in mm
brainSpace          = 'mni';                % native or mni
smoothnessDir       = [num2str(smoothKernelSize) 'mm-smoothed-' brainSpace 'space'];                     % Name of smoothed data directory

% if smoothing size is set to 0 use only realigned, slice time corrected
% and coregistered data 
dataDir = fullfile (derivesDir, softwareName, analysisPipeline, brainMask, glmUsed, smoothnessDir, 'SpecialMoment');  % data that is used for the analysis
destDir = fullfile (dataDir, 'Conjunction');  % where all the results of the TWO GLMs are stored
if ~isfolder(destDir)
    mkdir(destDir)
end

%% General settings
matlabbatch{1}.spm.util.imcalc.outdir = {destDir};
matlabbatch{1}.spm.util.imcalc.expression = '(i1 > 0) & (i2 > 0) & (i3 > 0)'; 
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm('defaults', 'FMRI');
%%%%%%%%%%%%%%%%%%%%%
% -log_10 p-values  %
% 3     => p=0.001  %
% 2.8   => p=0.005  %
% 2     => p=0.01   %
%%%%%%%%%%%%%%%%%%%%%
%% Conjunction 
% Effects > Controls
matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(dataDir, 'Appear vs Control Before/uncorrKclusterThrIMG.nii,1')
                                        fullfile(dataDir, 'Change vs Control Before/uncorrKclusterThrIMG.nii,1')
                                        fullfile(dataDir, 'Vanish vs Control Before/uncorrKclusterThrIMG.nii,1')
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'Conjunction_EffectsVsControl';
spm_jobman('run', matlabbatch);

% Pre vs Post
matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(dataDir, 'Appear Before vs Appear After/uncorrKclusterThrIMG.nii,1')
                                        fullfile(dataDir, 'Change Before vs Change After/uncorrKclusterThrIMG.nii,1')
                                        fullfile(dataDir, 'Vanish Before vs Vanish After/uncorrKclusterThrIMG.nii,1')
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'Conjunction_EffectsPreVsPost';
spm_jobman('run', matlabbatch);
% MagicvsControl - PrevsPost interaction
matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(dataDir, 'AppPre-ConPre vs AppPost-ConPost/uncorrKclusterThrIMG.nii,1')
                                        fullfile(dataDir, 'ChaPre-ConPre vs ChaPost-ConPost/uncorrKclusterThrIMG.nii,1')
                                        fullfile(dataDir, 'VanPre-ConPre vs Vanpost-ConPost/uncorrKclusterThrIMG.nii,1')
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'Conjunction_EffectsVsControl-PrevsPost_Interaction';
spm_jobman('run', matlabbatch);
