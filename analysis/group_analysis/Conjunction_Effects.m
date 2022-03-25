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
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "rawdata")\n\n'])
rootDir    = '/Users/vpl/Documents/Master_Thesis/DATA/MRI'; % uigetdir(homedir, 'Select Project Folder');
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
softwareName        = 'spm12';              % software used to create preprocessed data

% specify the name of the analysis pipeline
analysisPipeline    = 'spm12-sla';          % how is the folder named that contains first level results
brainMask           = 'WholeBrain';         % whole brain or ROI
glmUsed             = 'MagicEffects';
smoothKernelSize    = 9;                    % in mm
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
matlabbatch{1}.spm.util.imcalc.expression = '(i1 > 2.8) & (i2 > 2.8) & (i3 > 2.8)'; 
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm('defaults', 'FMRI');
%%%%%%%%%%%%%%%%%%%%%%%%%
% p-values              %
% T = 3.485 => p=0.001  %
% T = 2.8   => p=0.005  %
% T = 2.5   => p=0.01   %
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjunction 
% Effects > Controls
matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(dataDir, 'Appear > Control Before/spmT_0001.nii,1')
                                        fullfile(dataDir, 'Change > Control Before/spmT_0001.nii,1')
                                        fullfile(dataDir, 'Vanish > Control Before/spmT_0001.nii,1')
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'Conjunction_Effects>Control';
spm_jobman('run', matlabbatch);

% Pre > Post
matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(dataDir, 'Appear Before > Appear After/spmT_0001.nii,1')
                                        fullfile(dataDir, 'Change Before > Change After/spmT_0001.nii,1')
                                        fullfile(dataDir, 'Vanish Before > Vanish After/spmT_0001.nii,1')
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'Conjunction_EffectsPre>Post';
spm_jobman('run', matlabbatch);
% Magic>Control - Pre>Post interaction
matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(dataDir, 'AppPre-ConPre vs AppPost-ConPost/spmT_0001.nii,1')
                                        fullfile(dataDir, 'ChaPre-ConPre vs ChaPost-ConPost/spmT_0001.nii,1')
                                        fullfile(dataDir, 'VanPre-ConPre vs Vanpost-ConPost/spmT_0001.nii,1')
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'Conjunction_Effects>Control-Pre>Post_Interaction';
spm_jobman('run', matlabbatch);
