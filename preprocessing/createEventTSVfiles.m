
%% HEADER
% with the condition in a tsv file (actually in a txt file and then rename
% it to a tsv file) in BIDS conform order.
% The subjects are located in
% Master_Thesis/DATA/MRI/derivatives/PsychoPhysic all with the prefix sMag.
% For each run there are Block_XRun_X_<date-time>_{design, log, ptb} files
% containing all information about the experiment, conditions ratings,
% timings etc.
% log and ptb files are given to 'calculateHits.m' to get the ratings 

%% START
clear all; % first clear all variables, so we don't have any intervening variables.
tic; % start script.% Extract the surprise ratings from all subjects, save them in together

%% MISSING VALUES
% Unfortunately in some of our subjects the button for surprise=4 did not
% work. Judged by the notes taken by both experimenters in the following
% subjects a missing value is filled with a 4. 
% In cases where experimenters wrote something like "subject was sleepy,
% unconcentrated, bored etc" we do not change missing values

% TO-DO: discuss 
% -sub-10 
% -sub-12 --> in notes steht in run 11 hat einmal 5 gedrückt meinte aber 1
% korrigieren?
% -sub-14 notes sagen war müde (kaum n/a) und 4 in manchen runs detected
% --> würde sagen nein 
% -sub-18 --> gleiche wie bei sub-14
subsFillMisses = [1 8 9 11 12 15 17]; 



%% Define important details of your file structure and location
homedir = '/home/vplikat';
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "rawdata")\n\n'])
rootDir    = uigetdir(homedir, 'Select Project Folder');
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

psychoPhysicsDir = fullfile(derivesDir,'PsychoPhysic');
task = '_task-magic_';
formatSpec  = '%02i';

subPrefix         = 'sub-'; % input (['Please specify the prefix of your participant data in your SOURCE DATA.\n' ...
%'(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');
DICOMsubDirs       = spm_select('FPList', psychoPhysicsDir, 'dir', ['^' subPrefix]);  % ^ is needed so that the filter only searches for folders STARTING with DICOMprefix
DICOMsubDirs       = cellstr(DICOMsubDirs);  

for sub = 7%1:length(DICOMsubDirs)
    % Go into the PsychoPhysic directory of the current subject and get all
    % run data
    logs = cellstr(spm_select('FPList', DICOMsubDirs{sub}, 'log.mat'));
    ptbs = cellstr(spm_select('FPList', DICOMsubDirs{sub}, 'ptb.mat'));
    if length(logs) ~= length(ptbs)
        error('Not the same number of log files as ptb files')
    end
    
    % check if missing values have to be filled with value 4
    if any(ismember(subsFillMisses,sub))
        fillMissingValues = true;
    else
        fillMissingValues = false;
    end
    
    for r = 1:length(logs)
        % load the ptb and log file from the current run and give it to
        % 'calculateHits'. The returned Results are saved in a table
        % together with the corresponding condition
        % IMPORTANT: make sure emtpy values are saved as n/a (as it is
        % suggested by BIDS)
        load(logs{r},'log');
        load(ptbs{r},'ptb');
        [onset, duration, trialType, value, ratingOnset, ratingDuration] = calculateHits(log, ptb, fillMissingValues);
        event_struct.onset              = onset';
        event_struct.duration           = duration';
        event_struct.rating_onset       = ratingOnset';
        event_struct.rating_duration    = ratingDuration';
        event_struct.trial_type         = trialType';
        event_struct.value              = value';
%         event_struct.rating         = ratings';
%         event_struct.condition      = log.data.Condition;
%         event_struct.video_start    = log.data.VideoStart - log.data.VideoStart(1);
%         event_struct.rating_start   = log.data.Rating_vbl - log.data.VideoStart(1);
        event_table = struct2table(event_struct);
        fileName = ['sub-' num2str(sub,formatSpec) task 'run-' num2str(r,formatSpec) '_events'];
        writetable(event_table,fileName,'Delimiter','\t')
        [success,message] = movefile([fileName '.txt'], fullfile(rawDir, ['sub-' num2str(sub,formatSpec)], 'func', [fileName '.tsv']));
        if ~success
            warning(message)
        end
    end
end