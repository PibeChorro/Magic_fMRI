%% clear and close EVERYTHING
clear all;
close all;

%% Define important details of your file structure and location
homedir = '/home/vplikat';
workingDir = pwd;
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "derivatives")\n\n'])
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
softwareName        = 'snpm13';  
% get the data from the first level analysis pipeline you want
pipelineName        = 'spm12-fla';
brainMask           = 'WholeBrain';         % whole brain or ROI
conditionsAnalyzed  = 'MagicEffects';
smoothKernelSize	= 6;   % in mm
smoothKernelSpace   = 'mni';
% combine above specifications for a well structured file hierarchy
smoothnessDir       = [num2str(smoothKernelSize) 'mm-smoothed-' smoothKernelSpace 'space'];
flaDir              = fullfile(derivesDir, 'spm12', pipelineName, brainMask, conditionsAnalyzed, smoothnessDir);
wholeVideoGLMName   = 'WholeVideo';
specialMomentGLMName= 'SpecialMoment';
% specify the name of the processing pipeline
analysisPipeline    = 'snpm13-sla';

%% create a folder that contains the results of the second level analysis

secLevelDir = fullfile(derivesDir, softwareName, analysisPipeline, brainMask, conditionsAnalyzed, smoothnessDir);
if ~isfolder(secLevelDir)
    mkdir(secLevelDir)
end

%% Define what to do
do.specify      = 1;
do.compute      = 1;
do.inference    = 1;
% Which model to do
do.wholeVideo       = 0;
do.specialMoment    = 1;
% Which inference and correction
do.clusterInference = 1;
do.FWECorrection    = 0;

%% Settings
% specify format for folder numeration
formatSpec = '%04i';

folders = dir(fullfile(flaDir, specialMomentGLMName,'sub-*'));
subNames = {folders(:).name}; 

% load in one SPM.mat file to read out the contrasts
load(fullfile(flaDir, specialMomentGLMName,subNames{1},'SPM.mat'));
nContrasts = length(SPM.xCon);



%% The actual second level analysis
% Iterate over all contrasts
for C = 1:nContrasts 
    if do.specify
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignName = 'MultiSub: One Sample T test on diffs/contrasts';
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignFile = 'snpm_bch_ui_OneSampT';
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov = struct('c', {}, 'cname', {});
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.nPerm = 5000;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.vFWHM = [0 0 0];
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.ST.ST_later = -1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.im = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.em = {''};
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalc.g_omit = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 1;
        
        if do.wholeVideo
            dataDir         = fullfile(flaDir, wholeVideoGLMName);
            contrastsDirs   = {};

            for s = 1:length(subNames)

                currentDir              = fullfile(dataDir,subNames{s});
                contrastsDirs{end+1}    = fullfile(currentDir,['con_' num2str(C,formatSpec) '.nii']);

            end

            matlabbatch{1}.spm.tools.snpm.des.OneSampT.P = contrastsDirs';
            matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir = {fullfile(secLevelDir,specialMomentGLMName,SPM.xCon(C).name)};

            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);

        end

        %% special moment analysis
        if do.specialMoment
            dataDir         = fullfile(flaDir, specialMomentGLMName);
            contrastsDirs   = {};

            for s = 1:length(subNames)

                currentDir = fullfile(dataDir,subNames{s});
                contrastsDirs{end+1} = fullfile(currentDir,['con_' num2str(C,formatSpec) '.nii']);

            end

            matlabbatch{1}.spm.tools.snpm.des.OneSampT.P = contrastsDirs';
            matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir = {fullfile(secLevelDir,specialMomentGLMName,SPM.xCon(C).name)};

            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        clear matlabbatch;
    end

    if do.compute
        if do.wholeVideo
            matlabbatch{1}.spm.tools.snpm.cp.snpmcfg = {fullfile(secLevelDir,wholeVideoGLMName,SPM.xCon(C).name, 'SnPMcfg.mat')};
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end

        if do.specialMoment
            matlabbatch{1}.spm.tools.snpm.cp.snpmcfg = {fullfile(secLevelDir,specialMomentGLMName,SPM.xCon(C).name, 'SnPMcfg.mat')};
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        clear matlabbatch;
    end

    if do.inference
        if do.clusterInference
            %% cluster-level inference
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 0.001; %uncorrected p-value to define clusters
            if do.FWECorrection
                % FWE corrected
                matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05; % removes all clusters smaller than a threshold that will control the FWE rate
                matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'FWE_corr_clusterIMG';
%             elseif
%                 % uncorrected nonparametric p-value
%                 % WARNING1: produces p-values only for an a priori selected cluster
%                 % WARNING2: p-values obtained with an assumption of stationarity which is
%                 % likely to be inappropriate for VBM,MEG,EEG
%                 matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.PthC = 0.001; % p-value threshold to obtain uncorrected cluster inference
            else 
                % uncorrected with cluster size threshold
                matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.Cth = 10; 
                matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'uncorrKclusterThrIMG';
            end
        else
            %% voxel-level inference
            if do.FWECorrection
                % FWE corrected
                matlabbatch{1}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FWEth = 0.05;
                matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'FWE_corr_voxIMG';
%             elseif
%                 % uncorrected T-values
%                 matlabbatch{1}.spm.tools.snpm.inference.Thr.Vox.VoxSig.TFth = '3.5';
            else
                % uncorrected nonparametric p-value
                matlabbatch{1}.spm.tools.snpm.inference.Thr.Vox.VoxSig.Pth = 0.001;
                matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'uncorrVoxIMG';
            end
        end
        %% the rest
        matlabbatch{1}.spm.tools.snpm.inference.Tsign = 1;
        matlabbatch{1}.spm.tools.snpm.inference.Report = 'MIPtable';

        if do.wholeVideo
            matlabbatch{1}.spm.tools.snpm.inference.SnPMmat = {fullfile(secLevelDir,wholeVideoGLMName,SPM.xCon(C).name, 'SnPM.mat')};
            spm_jobman('run', matlabbatch);
        end

        if do.specialMoment
            cd (fullfile(secLevelDir,specialMomentGLMName,SPM.xCon(C).name))
            matlabbatch{1}.spm.tools.snpm.inference.SnPMmat = {fullfile(secLevelDir,specialMomentGLMName,SPM.xCon(C).name, 'SnPM.mat')};
            spm_jobman('run', matlabbatch);
        end
        clear matlabbatch;
    end
    %% save the created figure
    figHandle = findobj(allchild(0), 'flat', 'Type', 'figure');
    figName = 'Glasbrain.png';
    saveas(figHandle, fullfile(secLevelDir,specialMomentGLMName,SPM.xCon(C).name, figName), 'png')
end

cd (workingDir)
clear