% Supervisor: Pablo Grassi
% The following study is an fMRI experiment investigating the neural activation
% during the presentation of magic tricks.
% The experiment has a 3 (object) x 3 (trick) x 3 (effect) factor design
% Objects: ball, playingcard and stick
% Tricks: appear, vanish and colorchange
% Effect: magic (experime0ntal condition), nothing happens (first control)
% and surprise (second control)
% The experimental condition shows videos of magic tricks, the first
% control condition shows the same performance without the magical effect
% and the second control condition shows a similar performance with an
% unexpected event.
% The first controll condition is necessary to get similar visual input.
% The second controll condition serves to differentiate the neural
% correlate of supprise and violation of expectation.
% The experimental design is as followed:
% The experiment has three blocks. Each block is associated with one
% object. The blocks are divided into four runs. After the second run the
% tricks regarding the objects are revealed.3
% In each run the two control conditions are shown twice and the magic
% condition is shown four times.
% The order of presentation is randomized in every run and the order of
% blocks is randomized over the participants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MagicExperiment ()

try
    %% PREPARE EXPERIMENT
    % Settings for the Psychtoolbox
    ptb = MagicPTBSettings;
    
    %If true --> skips ET and fMRI triggers.
    log.ETdummy   = false;
    log.fMRIdummy = false;
    
    % Safety net for experimentator input. Only numbers (and test for sub) are
    % allowed
    
    checkinput = true;
    safetynet = true;
    while safetynet
        %% GET SUBJECT DATA
        log.subjNr = input('Enter subject Nr: ','s');
        log.block  = input('Which block [1,2,3]: ','s');
        log.run    = input('Which run [1-5]: ', 's');
        while checkinput
            
            if (str2double(log.subjNr)<=36 || strcmp (log.subjNr, 'test')) && str2double(log.block)<4 ...
                    && str2double (log.run)<6 && isnumeric(str2double(log.subjNr)) && ...
                    isnumeric(str2double(log.block)) && isnumeric(str2double(log.run))
                checkinput = false;
            else
                disp ('Invalid input. Make it right')
                log.subjNr = input('Enter subject Nr: ','s');
                log.block  = input('Which block [1,2,3]: ','s');
                log.run    = input('Which run [1-5]: ', 's');
            end
        end
        checkinput = true;
        while checkinput
            log.isGerman = input ('Does the subject understand german? ','s');
            if strcmp (log.isGerman, KbName(ptb.Keys.yes)) 
                checkinput = false;
            elseif strcmp (log.isGerman, KbName(ptb.Keys.no))
                checkinput = false;
            end
        end
        
        disp (['You specified subj ' log.subjNr ' in block ' log.block ' and run ' log.run])
        answer = input ('Are you sure you gave the right parameters?','s');
        if strcmp (answer, KbName(ptb.Keys.yes))
            safetynet = false;
        end
        checkinput = true;
    end
    %% Create unique filename
    c = clock;
    log.fileName = ['Data/Sub_' log.subjNr '/Block_' log.block 'Run_' log.run '_' num2str(c(2)) num2str(c(3)) num2str(c(4)) num2str(c(5))];
    
    if ~isfolder (['Data/Sub_' log.subjNr])
        mkdir (['Data/Sub_' log.subjNr])
    end
    %% Get design settings (variables, design, etc)
    % if it is a test, change subjNr to 1
    if strcmp(log.subjNr, 'test')
        log.subjNr ='1';
    end
    design = MagicDesignSettings(log);
    
    %% Preallocate Log structure
    log.data.Condition          = strings(length(design.Condition),1);
    log.data.Rating_vbl         = zeros(length(design.Condition),1);
    log.data.Rating_stimOn      = zeros(length(design.Condition),1);
    log.data.Rating_times       = zeros(length(design.Condition),1);
    log.data.Rating_missed      = zeros(length(design.Condition),1);
    log.data.Question_vbl       = zeros(length(design.Condition),1);
    log.data.Question_stimOn    = zeros(length(design.Condition),1);
    log.data.Question_times     = zeros(length(design.Condition),1);
    log.data.Question_missed    = zeros(length(design.Condition),1);
    log.data.VideoStart         = zeros(length(design.Condition),1);
    log.data.VideoEnd           = zeros(length(design.Condition),1);
    log.data.idDown   = [];
    log.data.timeDown = [];
    log.data.idUp   = [];
    log.data.timeUp = [];
    
    %% Preallocate Videos
    %For all runs apart from 3
    if log.run ~= '3'
        [design, tmpMoviePre]=preallocateVideos(design,ptb);
    end
    
    %% START EYETRACKER
    %Include Eyetracker
    if ~log.ETdummy
        [ET,log] = IncludeEyetracker(ptb,log);
        %% Calibration and Validation
        EyelinkDoTrackerSetup(ET);
        EyelinkDoDriftCorrection(ET);
        
        % start recording eye position
        Eyelink('StartRecording');
        % record a few samples before we actually start displaying
        WaitSecs(0.1);
        % mark zero-plot time in data file
        Eyelink('Message', 'SYNCTIME');
        % sets EyeUsed to -1 because we do not know which eye we use. This is
        % asked before we sample
        EyeUsed = -1;
    end
    
    %% AVOID KEYBOARD INPUT INTO MATLAB
    if ~ptb.Keys.debug
        ListenChar(-1); % Avoid input into Matlab windows only if you are not debugging.
    end
    
    %% SHOW INSTRUCTIONS
    if str2double(log.block)==1 && str2double(log.run) ==1
        %...........................Show instructions in very first run.............................%
        DrawFormattedText (ptb.window, design.Introduction, 'center', 'center',ptb.FontColor);
        Screen ('Flip', ptb.window);
        KbWait(ptb.Keys.kbrd2);
        
        DrawFormattedText (ptb.window,design.Instruction1, 'center', 'center', ptb.FontColor);
        Screen ('Flip', ptb.window);
        WaitSecs (0.5); KbWait;
        
        DrawFormattedText (ptb.window,design.Instruction2, 'center', 'center', ptb.FontColor);
        Screen ('Flip', ptb.window);
        WaitSecs (0.5); KbWait;
        
        DrawFormattedText (ptb.window,design.Instruction3, 'center', 'center', ptb.FontColor);
        Screen ('Flip', ptb.window);
        WaitSecs (0.5);KbWait;
        
    elseif str2double(log.run) == 3 % Show revealing instructions
        DrawFormattedText (ptb.window, design.RevealingsWillBeShown, 'center', 'center', ptb.FontColor);
        Screen ('Flip', ptb.window);
        WaitSecs (0.5);KbWait;
    else % Only the last one
        DrawFormattedText (ptb.window,design.Instruction3, 'center', 'center', ptb.FontColor);
        Screen ('Flip', ptb.window);
        WaitSecs (0.5);KbWait;
        %............................................................................................%
    end
    
    %% Stop and remove events in queue
    KbQueueStop(ptb.Keys.kbrd2);
    KbEventFlush(ptb.Keys.kbrd2);
    KbQueueStop(ptb.Keys.kbrd1);
    KbEventFlush(ptb.Keys.kbrd1);
    
    %% START FMRI TRIGGERS
    %Include fMRIScanner
    if ~log.fMRIdummy && log.run ~= '3'
        [ptb, log] = MagicGetTrigger(ptb,log,design.nDummies);
    end
    
    %% Start KbQueues
    KbQueueStart(ptb.Keys.kbrd2); % Subjects
    KbQueueStart(ptb.Keys.kbrd1); % Experimentors
    
    %% SHOW VIDEOS AND START MEASURING!
    if  log.run ~= '3'
        for video = 1:length(design.Condition) % Number of Videos per run
            tic
            if ~log.ETdummy
                trialNumer = [log.block log.run video];
                Eyelink('Message', sprintf('TRIAL %u', trialNumer));
            end
            
            %.....SHOW VIDEO....%
            %fprintf ('Presenting Video %s \n', design.Condition(video)); % display the video shown
            log.data.VideoStart(video) = GetSecs(); % get start time
            ShowVideo(ptb, tmpMoviePre(video), design.Condition(video), log.run); % show video
            log.data.VideoEnd(video)   = GetSecs(); % get end time
            %...................%
            
            %...GET RESPONSES...%
            % In the ultimate FINAL version we do not ask if they have seen
            % a magic trick or not
%             DrawFormattedText (ptb.window,design.WasMagic,'center','center',ptb.FontColor);
%             [log.data.Question_vbl(video), log.data.Question_stimOn(video), log.data.Question_times(video),...
%                 log.data.Question_missed(video)]  = Screen ('Flip', ptb.window, ...
%                 log.data.VideoStart(video) + ptb.VideoLength);
            
            DrawFormattedText (ptb.window,design.HowSurprising,'center','center',ptb.FontColor);
            [log.data.Rating_vbl(video), log.data.Rating_stimOn(video), log.data.Rating_times(video), ...
                log.data.Rating_missed(video)] = Screen ('Flip', ptb.window, ...
                log.data.VideoStart(video) + ptb.VideoLength);%log.data.Question_vbl(video) + ptb.ResponseTime);
            
            WaitSecs(ptb.ResponseTime); % wait for answers
            
            %.....Write data in .mat file......%
            log.data.Condition(video)       = design.Condition(video);
            toc
            %... Check if Experimentors pressed escape ...%
            [~, firstPress] = KbQueueCheck(ptb.Keys.kbrd1);
            if firstPress(ptb.Keys.escape)
                log.end = 'Escape'; % Finished by escape key
                log.ExpNotes = input('Notes:','s');
                savedata(log,design,ptb);
                cleanup(log);
                warning('Experiment was terminated by Escapekey');
                return;
            end
        end
        
        %...............Show final text...........%
        if log.run ~= '5'
            DrawFormattedText (ptb.window,design.RunIsOver, 'center', 'center', ptb.FontColor);
        else
            DrawFormattedText (ptb.window,design.BlockIsOver, 'center', 'center', ptb.FontColor);
        end
        Screen ('Flip', ptb.window);
        WaitSecs (1); KbWait;
        %............................REVEALING....................................%
    elseif log.run == '3'
        
        answer = 0;
        tmpTrickFile = strcat([pwd,'/../Stimuli/'], design.Tricks(:),'.mp4');
        tmpRevealingFile = strcat([pwd,'/../Stimuli/'], design.Revealing(:),'.mp4');
        for i = 1:length (tmpTrickFile)
            TrickFiles{i} = char(tmpTrickFile(i,:));
            RevealingFiles{i} = char(tmpRevealingFile(i,:));
        end
        
        log.reveals = zeros(length(design.Revealing),1); % Reveals counter
        
        %..............Show all revealing videos..................%
        for r = 1:length (design.Revealing)
            
            answer = 0;
            while answer~= 1
                fprintf ('Presenting Revealing %s\n', design.Revealing(r));
                %..........................VIDEOS.............................%
                % Add a reveal to counter
                log.reveals(r) = log.reveals(r)+1;
                
                % Show trick
                DrawFormattedText (ptb.window, ['Trick: ' num2str(r)], 'center', 'center');
                Screen ('Flip', ptb.window); WaitSecs (0.5);
                
                tmpTrickPre = Screen ('OpenMovie', ptb.window, TrickFiles{r});
                ShowVideo(ptb,tmpTrickPre,design.Tricks(r), log.run);
                
                % Show revealing video
                DrawFormattedText (ptb.window, ['Aufloesung: ' num2str(r)], 'center', 'center');
                Screen ('Flip', ptb.window);  WaitSecs (0.5);
                
                tmpRevealingPre = Screen ('OpenMovie', ptb.window, RevealingFiles{r});
                ShowVideo(ptb,tmpRevealingPre,design.Revealing(r), log.run);
                
                %........................RESPONSES............................%
                % Get responses. If yes, continue to next video, if no, repeat.
                % Answer IS needed.
                DrawFormattedText (ptb.window, design.TrickUnderstood, 'center', 'center');
                Screen ('Flip', ptb.window);
                
                while true
                    [KeyIsPressed, presstime, KeyCode] = KbCheck(ptb.Keys.kbrd2);
                    if KeyIsPressed
                        if find(KeyCode) == ptb.Keys.yes_sub
                            answer = 1;
                            break
                        elseif find(KeyCode) == ptb.Keys.no_sub
                            answer = 0;
                            break
                        end
                    end
                end
                
                % Cancel experiment if experimentors pressed ESCAPE
                [~, firstPress] = KbQueueCheck(ptb.Keys.kbrd1);
                if firstPress(ptb.Keys.escape)
                    log.end = 'Escape'; % Finished by escape key
                    log.ExpNotes = input('Notes:','s');
                    savedata(log,design,ptb);
                    cleanup(log);
                    warning('Experiment was terminated by Escapekey');
                    return;
                end
            end
        end
        
        %...............Show final text...........%
        DrawFormattedText (ptb.window,design.RunIsOver, 'center', 'center', ptb.FontColor);
        Screen ('Flip', ptb.window);
        WaitSecs (1); 
    end
    
    %................Save files...............%
    log.end = 'Success'; % Add a last summary.
%     log.ExpNotes = input('Notes:','s');
    savedata(log,design,ptb);
    cleanup(log);
    
catch MY_ERROR
    log.end = 'Finished with errors';
    log.error = MY_ERROR;
    log.ExpNotes = input('Notes:','s');
    savedata(log,design,ptb);
    cleanup(log);
    rethrow (MY_ERROR);
end
end

function [design,tmpMoviePre] =preallocateVideos(design,ptb)

tmpMovieFile = strcat([pwd,'/../Stimuli/'], design.Condition(:),'.mp4');
for i = 1:length (tmpMovieFile)
    MovieFiles{i} = char(tmpMovieFile(i,:));
    tmpMoviePre(i) = Screen ('OpenMovie', ptb.window, MovieFiles{i});
end

clear tmpCopy % just as a sanity check
for i = 1:length(design.Condition)
    tmpCopy(i) = design.Condition(i); 
    if ~contains(tmpCopy(i),'Control') 
        if contains(tmpCopy(i),'Vanish1_Magic') && sum(contains(tmpCopy,'Vanish1_Magic'))==1 
            tmpCopy(i) = strcat (tmpCopy(i), "_F");
        elseif contains(tmpCopy(i),'Vanish2_Magic') && sum(contains(tmpCopy,'Vanish2_Magic'))==2
            tmpCopy(i) = strcat (tmpCopy(i), "_F");
        elseif contains(tmpCopy(i),'Appear1_Magic') && sum(contains(tmpCopy,'Appear1_Magic'))==1
            tmpCopy(i) = strcat (tmpCopy(i), "_F");
        elseif contains(tmpCopy(i),'Appear2_Magic') && sum(contains(tmpCopy,'Appear2_Magic'))==2
            tmpCopy(i) = strcat (tmpCopy(i), "_F");
        elseif contains(tmpCopy(i),'Change1_Magic') && sum(contains(tmpCopy,'Change1_Magic'))==1
            tmpCopy(i) = strcat (tmpCopy(i), "_F");
        elseif contains(tmpCopy(i),'Change2_Magic') && sum(contains(tmpCopy,'Change2_Magic'))==2
            tmpCopy(i) = strcat (tmpCopy(i), "_F");
        elseif logical(contains(tmpCopy(i),'Surprise1')) || logical(contains(tmpCopy(i),'Surprise3'))
            tmpCopy(i) = strcat (tmpCopy(i), "_F");
        end
    end  
end

tmpCopy = reshape (tmpCopy, [length(tmpCopy),1]);
design.Condition = tmpCopy;
end
