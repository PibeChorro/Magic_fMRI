%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for the first level analysis of the "Magic" Experiment.
% Idea for study: Dr. Pablo R Grassi in cooperation with Prof. Dr. Andreas
% Bartels in the 'Vision and Cognition Lab' at the University of Tübingen.
% Development, implementation and execution: Vincent Plikat and Dr. Pablo R
% Grassi in cooperation with Prof. Dr. Adreas Bartels
% Started August 2019 as the Master Thesis of Vincent Plikat
% Data: functional MRI, pupil dilation, gaze position and behavioural
% rating (1-5). In this analysis only the fMRI data will be analyzed.
% Experimental design: The Experiment was divided into 3 blocks. Each block
% consisted of 4 experimental runs. In each run subjects viewed 24 videos
% (each video is considered a trial) of three different categories: Magic
% Control and Surprise.
% The videos in each block were associated with one object (Balls, Cards
% and Sticks)
% There were 6 Magic videos with 6 coresponding control videos and 3
% surpise videos per block. 3 Different magic effects were used: Object
% appearing, vanishing and color change. The magic and surprise videos were
% presented twice, the control videos only once. (2*6 magic + 6 control + 2*3
% surpsise = 24 videos). Half of all video presentations were flipped along
% the y-axis and the order of videos was pseudorandomized in a way that the
% same video was never presented in two consecutive trials.
% The videos were exactly 14 seconds long and followed by a 2 second
% answering phase in which the subject was asked to rate how surprising the
% video's content was.
% After the second run in each block the underlying methods behind each
% magic trick was presented (during this period we did not measure fMR
% data). Only when every trick was understood the following two runs were
% presented
% --> a 2(pre post revelation)*3(objects)*3(magic effect)*3(video type)
% design.
% In this analysis we analyze the data of the whole video presentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; % first clear all variables, so we don't have any intervening variables.
tic; % start script.

%% Define important details of your file structure and location
homedir = '/home/vplikat';
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "rawdata")\n\n'])
rootDir    =  uigetdir(homedir, 'Select Project Folder');
if rootDir == 0
    error('No folder was selected --> I terminate the script')
end

% Set sourcedata directory. This is needed to get DICOM header information
% and the names of subjects in the PsychoPhysics data folder
sourceDir  = fullfile(rootDir,'sourcedata');
if ~isfolder(sourceDir)
    fprintf(['It appears you do not have a "sourcedata" folder.\n'...
        'Please select the folder that contains your DICOMS.'])
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
analysisPipeline    = 'spm12-fla';          % how is the folder named that contains first level results
brainMask           = 'WholeBrain_24MP';         % whole brain or ROI
conditionsAnalyzed  = 'VideoTypes';         % Magic, Control and Surprise videos are one regressor each (Response as well, but that is less important)
smoothKernelSize	= 6;                    % in mm
smoothKernelSpace   = 'mni';                % mni or native (mni makes more sense, native rather for explorative analysis... maybe)
% combine above specifications for a well structured file hierarchy
smoothnessDir       = [num2str(smoothKernelSize) 'mm-smoothed-' smoothKernelSpace 'space'];                     % Name of smoothed data directory
% Get name, location and number of sourcedata subjects

DICOMprefix         = 'sMag'; % input (['Please specify the prefix of your participant data in your SOURCE DATA.\n' ...
%'(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');
DICOMsubNames       = spm_select('List', sourceDir, 'dir', ['^' DICOMprefix]);  % ^ is needed so that the filter only searches for folders STARTING with DICOMprefix
DICOMsubNames       = cellstr(DICOMsubNames);                                   % cellstring format is needed for spm

% Get name, location and number of preprocessed subjects
pipelineName        = 'spm12-preproc';                                                  % how is the folder named that contains preprocessed data
realignedDir        = fullfile(derivesDir, softwareName, pipelineName, 'realigned');    % needed for realignment files as regressors of no interest
if smoothKernelSize == 0
    destDir = fullfile (derivesDir, softwareName, analysisPipeline, brainMask, conditionsAnalyzed, [smoothKernelSpace 'space']);  % where all the results of the TWO GLMs are stored
    if strcmp(smoothKernelSpace,'native')
        dataDir = fullfile (derivesDir, softwareName, pipelineName, 'coregistered');  % data that is used for the analysis
    elseif strcmp(smoothKernelSpace,'mni')
        dataDir = fullfile (derivesDir, softwareName, pipelineName, 'normalized');  % data that is used for the analysis
    end
else
    destDir = fullfile (derivesDir, softwareName, analysisPipeline, brainMask, conditionsAnalyzed, smoothnessDir);  % where all the results of the TWO GLMs are stored
    dataDir = fullfile (derivesDir, softwareName, pipelineName, smoothnessDir);  % data that is used for the analysis
end
subNames            = spm_select('List', dataDir, 'dir', 'sub-');
subNames            = cellstr(subNames);

%% Multiband Factor
multibandFactor = 2;

%% Get the information about the videos
% read in the .mat file, that contains the information about the magic
% moment
fprintf('Please select the .mat file, that contains the information about the special moments.\n\n')
videoInfoMatfileDir = uigetfile(pwd, 'Select the .mat file');
if videoInfoMatfileDir == 0
        error('No .mat file was selected --> I terminate the script')
end
load(videoInfoMatfileDir);

% information about the videos. Important is only the framerate and thus
% the time every frame was presented
fps         = 25;
frameTime   = 1/fps;

%% Define what to do
do.SpecifyDesign    = 1;
do.loadlog          = 1; % load LOG files!
do.estimate         = 1;
do.DefContrasts     = 1;
% Which model to do
do.wholeVideo       = 1;
do.specialMoment    = 1;


%% Get experimental design parameters
fla.conditionNames  = {
    'Magic';    'Control';...
    'Surprise'; 'Response'...
    };
magicEffects = {
    'Appear1'   ; 'Appear2';    ...
    'Vanish1'   ; 'Vanish2';    ...
    'Change1'   ; 'Change2'     ...
    };

% how many conditions we have (in total and per "type")
fla.numConditions   = length(fla.conditionNames);
numMag              = 6;
numCon              = 6;
numSur              = 1;
numRaPara           = 24; % default number of realignment parameter, in model estimation overwritten but still 24
% how many blocks we had
numBlocks           = 3;
% should the movement be used as regressors of no interest
fla.realignmentParametersFlag  = 1;

for s = 1:length(subNames)
    %% Define where to look for functional MRI data and the logs that contain information about stimulus on/offsets
    smoothedDataDir     = fullfile(dataDir,         subNames{s}, 'func');
    realignedDataDir    = fullfile(realignedDir,    subNames{s}, 'func');
    psyphysicDataDir    = fullfile(derivesDir, 'PsychoPhysic',subNames{s});
    % Further information - number of runs and where a DICOM file can be
    % found
    runs                = cellstr(spm_select('List', smoothedDataDir, '.nii')); 
    numRuns             = 12;%length(runs); 
    sourcedataRuns      = spm_select('FPList',fullfile(sourceDir,DICOMsubNames{s},'func'), 'dir','^run*');
    sourcedataRuns      = cellstr(sourcedataRuns);
    %% load a dicom header that contains information needed for analysis
    dicomFiles  = spm_select('FPList', sourcedataRuns{1}, '**.IMA');    % ** is needed because one * stops at the first period. ** filters for the whole filename
    dicomFiles  = cellstr(dicomFiles);
    hdr         = spm_dicom_headers(dicomFiles{1});
    
    %% Model specification of FLA
    if do.SpecifyDesign == 1 % Model specification
        
        %% DEFINE MODEL PARAMETERS GENERAL
        matlabbatch{1}.spm.stats.fmri_spec.timing.units     = 'secs';                       % 'secs' or 'scans' unit
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT        = hdr{1}.RepetitionTime/1000;   % TR in seconds!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;70bf5c7f.0904
        % suggest to set the microtime resolution to the number of slices,
        % if slice time correction was implemented and microtime onset to
        % the half. 
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t    = hdr{1}.Private_0019_100a/multibandFactor;             % Number of slices - before it was 16
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0   = round(hdr{1}.Private_0019_100a/(multibandFactor*2));	% Number of slices - before it was 8 --> T0 = reference slice is the one in the middle. 18 or 8
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];                % derivatives
        matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;                    % CHECK
        matlabbatch{1}.spm.stats.fmri_spec.global           = '';                   % Global scaling? if so: 'Scaling' if not empty: ''
        matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;                  % Masking threshold
        matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};                 % Mask?
        matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
        
        %% Get the realigment files
        AlignmentFiles = spm_select('FPList',realignedDataDir,'^rp_.*.txt');
        
        % Load .mat files that contain condition information, such as
        % condition names, onsets, offsets, etc.
        if do.loadlog == 1 % get the mat files
            logMatfiles = spm_select('FPList', psyphysicDataDir, 'log.mat'); % SELECT MAT FILES
        end
        
        
        for r = 1:numRuns % For each run
            
            % Find functional files
            dirfiles    = spm_select('ExtFPList',smoothedDataDir, runs(r), Inf);
            
            % check if any files were selected. If not stop the script
            if strcmp(dirfiles,'')
                warning('No files selected!');
                return;
            else
                nfiles     = length(dirfiles);
                fprintf('We got %s files! \n', num2str(nfiles));
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans        = cellstr(dirfiles); %For every run scans are added to the matlabbatch
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi        = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg    = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf          = 128;        % High-pass filter
            
            % delete the variable 'dirfiles'
            clear dirfiles
            
            %......Include realignment parameters for each run
            if fla.realignmentParametersFlag == 1
                
                fprintf('Adding realignment parameters to design matrix! \n');
                
                % Define names for 'Motion Regressors' (aka. Realignment
                % Parameters)
                raParamNames = [
                    'rp01'; 'rp02'; 'rp03';...  % translation
                    'rp04'; 'rp05'; 'rp06';...  % rotation
                    'rp07'; 'rp08'; 'rp09';...  % translation 1st derivative
                    'rp10'; 'rp11'; 'rp12';...  % rotation 1st derivative
                    'rp13'; 'rp14'; 'rp15';...  % translation square
                    'rp16'; 'rp17'; 'rp18';...  % rotation square
                    'rp19'; 'rp20'; 'rp21';...  % translation 1st derivative sqaure
                    'rp22'; 'rp23'; 'rp24'      % rotation 1st derivative square
                    ];
                numRaPara = length(raParamNames);
                
                % Load alignment parameters
                raParamValues       = load(AlignmentFiles(r,:));    % get realignment parameters
                [~, raDeriValues]   = gradient(raParamValues);      % the first return value is the derivative along x axis
                nuis                = zscore([raParamValues raDeriValues]); 
                quadterms           = nuis.^ 2;
                all_raParams        = [nuis zscore(quadterms)];
                % Add realignment parameters as regressors of no interest.
                for P = 1:length (raParamNames)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(P).name = raParamNames(P,:);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(P).val  = all_raParams(:,P);
                end
            end
        end
        
        if do.wholeVideo
            %% Name, create and set directory for this analysis
            analysisName    = 'WholeVideo';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});

            if ~isfolder(betaLoc)
                mkdir(betaLoc);
            end

            matlabbatch{1}.spm.stats.fmri_spec.dir = {betaLoc};  % The directory, the SPM.mat and all betas are written

            for r = 1:numRuns % For each run

                if do.loadlog == 1

                    % IMPORTANT: this throws an error if there is a space in
                    % your path.
                    load(strtrim(logMatfiles(r,:)));

                    % specify the trial onsets and durations.
                    % We iterate over our conditions specified above
                    fla.trialOnset={}; fla.trialDuration={};
                    for reg = 1:length (fla.conditionNames) - 1 % The last condition is the response. We set it manually at the end
                        % from the matlab log file select the condition names.
                        % Look which index contains the current condition.
                        % Those indices are used to select the trial onset from
                        % the matlab log file. The very first EXPERIMENTAL
                        % trigger is substracted from Video start, to set
                        % the measurement onset = 0
                        fla.trialOnset{end+1}      = log.data.VideoStart(contains(log.data.Condition,fla.conditionNames(reg))) - log.Keys.trig(end);

                        % The same procedure as above to get the indices, but
                        % for rating onset (which is stim offset) minus stim
                        % onset
                        fla.trialDuration{end+1}   = log.data.Rating_stimOn(contains(log.data.Condition,fla.conditionNames(reg)))...
                            -log.data.VideoStart(contains(log.data.Condition,fla.conditionNames(reg)));
                    end
                    % We just add 2 seconds for the response in the last entry of trialDuration and the end of each Video as trialOnset
                    fla.trialOnset{end+1}         = log.data.Rating_stimOn - log.Keys.trig(end);
                    fla.trialDuration{end+1}      = ones(1, length (fla.trialOnset{end}))*2;
                else % don't get the matfile log files.

                    %PLACE HOLDER

                end

                % Condition names:
                for cc = 1:fla.numConditions % for CurrentCondition
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).name        = fla.conditionNames{cc};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).onset       = fla.trialOnset{cc};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).duration    = fla.trialDuration{cc};

                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).pmod        = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).tmod        = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).orth        = 1;
                end
            end
        
        %% Specify design in SPM.mat
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        
        if do.specialMoment
            %% Name, create and set directory for this analysis
            analysisName    = 'SpecialMoment';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});

            if ~isfolder(betaLoc)
                mkdir(betaLoc);
            end

            matlabbatch{1}.spm.stats.fmri_spec.dir = {betaLoc};   % The directory, the SPM.mat and all betas are written

            for r = 1:numRuns % For each run

                if do.loadlog == 1

                    % IMPORTANT: this throws an error if there is a space in
                    % your path. 
                    load(strtrim(logMatfiles(r,:)));

                    % specify the trial onsets and durations.
                    % We iterate over Magic and Control conditions 
                    % specified above
                    fla.trialOnset={}; fla.trialDuration={};
                    %% Magic and Control
                    for reg = 1:length (fla.conditionNames) - 2 % the last two conditions are Surprise and Rating
                        tmp = []; % An empty array to contain specialmoment onsets all videos belonging to one condition
                        for magEff = 1:length(magicEffects)
                            % Get the indices for for every video belonging
                            % to the current condition
                            VideoNames = log.data.Condition(contains(log.data.Condition,fla.conditionNames(reg)) & contains(log.data.Condition, magicEffects(magEff)));
                            % The videoname may contain an '_F' for the flip
                            % condition. Therefore we need to get the normal video
                            % name
                            if contains(VideoNames{1},'_F')
                                VideoName   =VideoNames{1};
                                VideoName   =VideoName(1:end-2);
                            else
                                VideoName   =VideoNames{1};
                            end
                            % Select the frame of the special moment from the
                            % 'do' struct
                            SpecialMoment               = do.both_frames_of_effect(contains(do.ListOfVideos,VideoName));
                            SpecialMomentOnset          = SpecialMoment*frameTime; % multiply with the frameduration

                            % The very first EXPERIMENTAL trigger is 
                            % substracted from Video start, to set the
                            % measurement onset = 0
                            tmp = [tmp; log.data.VideoStart(contains(log.data.Condition,VideoName))-log.Keys.trig(end) + SpecialMomentOnset];
                        end
                       
                        fla.trialOnset{end+1}       = tmp;
                        fla.trialDuration{end+1}    = zeros(1,length(tmp));
                    end
                    
                    %% Surprise
                    % Do the same as in the for loop above, but for the
                    % three surprise videos
                    reg = 3;
                    tmp = [];
                    
                    for surVer = 1:3
                        VideoNames = log.data.Condition(contains(log.data.Condition,fla.conditionNames(reg)) & contains(log.data.Condition, num2str(surVer)));
                        % The videoname may contain an '_F' for the flip
                        % condition. Therefore we need to get the normal video
                        % name
                        if contains(VideoNames{1},'_F')
                            VideoName   =VideoNames{1};
                            VideoName   =VideoName(1:end-2);
                        else
                            VideoName   =VideoNames{1};
                        end
                        % Select the frame of the special moment from the
                        % 'do' struct
                        SpecialMoment               = do.both_frames_of_effect(contains(do.ListOfVideos,VideoName));
                        SpecialMomentOnset          = SpecialMoment*frameTime; % multiply with the frameduration
                        
                        % The very first EXPERIMENTAL trigger is 
                        % substracted from Video start, to set the
                        % measurement onset = 0
                        tmp = [tmp; log.data.VideoStart(contains(log.data.Condition,VideoName))-log.Keys.trig(end) + SpecialMomentOnset];
                    end

                    fla.trialOnset{end+1}       = tmp;
                    fla.trialDuration{end+1}    = zeros(1,length(tmp));
                    %% Response
                    % We just add 2 seconds for the response in the last entry of trialDuration and the end of each Video as trialOnset
                    fla.trialOnset{end+1}       = log.data.Rating_stimOn - log.Keys.trig(end);
                    fla.trialDuration{end+1}    = ones(1, length (fla.trialOnset{end}))*2;
                else % don't get the matfile log files.

                  %PLACE HOLDER

                end % End load mat files

                % Condition names:
                for cc = 1:fla.numConditions % for CurrentCondition
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).name     = fla.conditionNames{cc};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).onset    = fla.trialOnset{cc};
                    % For this analysis we do  set the duration to 0, since we
                    % are only interested in the effect of the special onset
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).duration = fla.trialDuration{cc};

                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).pmod        = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).tmod        = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).orth        = 1;
                end
            end
        
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        clear matlabbatch;
        
    end % model specification
    
    %% Estimate design:
    if do.estimate == 1
        % General settings
        matlabbatch{1}.spm.stats.fmri_est.write_residuals  = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        if do.wholeVideo
            % Set directory for whole Video analysis
            analysisName    = 'WholeVideo';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});

            matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(betaLoc, 'SPM.mat')};
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        
        if do.specialMoment
            % Set directory for special moment analysis
            analysisName    = 'SpecialMoment';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});

            matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(betaLoc, 'SPM.mat')};
        
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        clear matlabbatch;
    end % end estimate
    
    
    %% Define contrast
    if do.DefContrasts == 1
        % General settings
        % Contrast Names:
        ContrastNames = ...
            {
            'Magic vs Control';... %1
            'Control vs Magic';... %2
            'Magic vs Control Before';               ... %3
            'Control vs Magic Before';               ... %4
            'Magic Before vs Magic After';           ... %5
            'Magic After vs Magic Before';           ... %6
            'Magic vs Surprise Before';              ... %7
            'Surpise vs Magic Before';               ... %8
            'Surprise vs Control';                   ... %9
            'Control vs Surprise';                   ... %10
            'Magic vs Control After';                ... %11
            'Control vs Magic After';                ... %12
            'Magic vs Surprise After';               ... %13
            'Surpise vs Magic After';                ... %14
            'MagPre-ConPre vs MagPost-ConPost';     ... %15
            'MagPost-ConPost vs MagPre-ConPre';     ... %16
            % Contrats to outrule the timeconfound by comparing run 1vs2
            % and run 2vs3 - the same time difference, but the first is pre
            % vs pre and the other is pre vs post revelation
            'Magic PreVsPre (run 1vs2)';            ... %17
            'Magic PreVsPost (run 2vs3)';           ... %18
            % Additional contrast to outrule the timeconfound by comparing
            % Control pre vs post
            'Control Before vs Control After';       ... %19
            'Video vs Response';                    ... %20
            'Response vs Video'                     ... %21
            };
        
        % Contrast values
        %       PreRevelation Magic Control Surprise    Response    Realignment             PostRevelation  Magic videos    Control Surprise    Response   Realignment
        C1 = repmat ([repmat([1     -1      0           0           zeros(1,numRaPara)], 1,2)       repmat([1               -1      0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C2 = repmat ([repmat([-1     1      0           0           zeros(1,numRaPara)], 1,2)       repmat([-1              1       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C3 = repmat ([repmat([1     -1      0           0           zeros(1,numRaPara)], 1,2)       repmat([0               0       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C4 = repmat ([repmat([-1     1      0           0           zeros(1,numRaPara)], 1,2)       repmat([0               0       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C5 = repmat ([repmat([1      0      0           0           zeros(1,numRaPara)], 1,2)       repmat([-1              0       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C6 = repmat ([repmat([-1     0      0           0           zeros(1,numRaPara)], 1,2)       repmat([1               0       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C7 = repmat ([repmat([1      0     -1           0           zeros(1,numRaPara)], 1,2)       repmat([0               0       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C8 = repmat ([repmat([-1     0      1           0           zeros(1,numRaPara)], 1,2)       repmat([0               0       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C9 = repmat ([repmat([0     -1      1           0           zeros(1,numRaPara)], 1,2)       repmat([0               -1      1           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C10= repmat ([repmat([0      1      -1          0           zeros(1,numRaPara)], 1,2)       repmat([0               1       -1          0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C11= repmat ([repmat([0      0      0           0           zeros(1,numRaPara)], 1,2)       repmat([1               -1      0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C12= repmat ([repmat([0      0      0           0           zeros(1,numRaPara)], 1,2)       repmat([-1              1       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C13= repmat ([repmat([0      0      0           0           zeros(1,numRaPara)], 1,2)       repmat([1               0       -1          0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C14= repmat ([repmat([0      0      0           0           zeros(1,numRaPara)], 1,2)       repmat([-1              0       1           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C15= repmat ([repmat([1     -1      0           0           zeros(1,numRaPara)], 1,2)       repmat([-1              1       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        C16= repmat ([repmat([-1     1      0           0           zeros(1,numRaPara)], 1,2)       repmat([1               -1      0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);
        %   FirstRun Magic  Control Surprise    Response    Realignment     SecondRun   Magic   Control Surprise    Response    Realignment Postrevelation  Magic   Control Surprise    Response    Realignment
        C17= repmat ([1     0       0           0           zeros(1,numRaPara)          -1      0       0           0           zeros(1,numRaPara) repmat([0       0       0           0           zeros(1,numRaPara)], 1,2)],1,numBlocks);
        %   FirstRun Magic  Control Surprise    Response    Realignment SecondRun   Magic   Control Surprise    Response    Realignment   ThirdRun  Magic   Control Surprise    Response    Realignment FourthRun   Magic   Control Surprise    Response    Realignment
        C18= repmat ([0     0       0           0           zeros(1,numRaPara)      1       0       0           0           zeros(1,numRaPara)      -1      0       0           0           zeros(1,numRaPara)      0       0       0           0       	zeros(1,numRaPara)],1,numBlocks);
        %       PreRevelation Magic videos  Control Surprise    Response    Realignment             PostRevelation  Magic videos    Control Surprise    Response   Realignment
        C19= repmat ([repmat([0              1      0           0           zeros(1,numRaPara)], 1,2)       repmat([0              -1       0           0           zeros(1,numRaPara)], 1,2)], 1, numBlocks);

        %                       Magic No Magic    Surprise    Response    Realigment
        C20= repmat (repmat([   1     1           1           -3          zeros(1,numRaPara)],1,4),1,numBlocks);
        C21= repmat (repmat([   -1    -1          -1           3          zeros(1,numRaPara)],1,4),1,numBlocks);
        % Combine all Contrasts in one Matrix
        Contrasts = [C1; C2; C3; C4; C5; C6; C7; C8; C9; C10; C11; C12; C13; C14; C15; C16; C17; C18; C19; C20; C21];
        
        % safety net: check if sum of contrasts is 0
        if any(sum(Contrasts,2))
            error('One of the contrasts sum is not zero')
        end
        
        %% Define contrasts:
        
        matlabbatch{1}.spm.stats.con.delete = 1; % Delete existing Contrasts
        
        % Write in spm-struct:
        for C = 1:size(Contrasts, 1)
            matlabbatch{1}.spm.stats.con.consess{C}.tcon.name       = ContrastNames{C,:};
            matlabbatch{1}.spm.stats.con.consess{C}.tcon.sessrep    = 'none'; %'replsc' if repeat for sessions
            matlabbatch{1}.spm.stats.con.consess{C}.tcon.convec     = Contrasts(C,:); % or weights?
        end
        
        if do.wholeVideo
            % Set directory for whole Video analysis
            analysisName    = 'WholeVideo';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});
            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(betaLoc, 'SPM.mat')};
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        
        %% Contrasts for analysis of the special moment. 
        
        if do.specialMoment
            % Set directory for special moment analysis
            analysisName    = 'SpecialMoment';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});
            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(betaLoc, 'SPM.mat')};
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        
        clear matlabbatch;
        
    end % define contrast
    
    toc;
    
end

clear all