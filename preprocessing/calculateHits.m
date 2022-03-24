function [onset, duration, trialType, value, ratingOnset, ratingDuration] = calculateHits(log,ptb, fillMisses)

onset           = [];
duration        = [];
trialType       = {};
value           = {};
ratingOnset     = [];
ratingDuration  = [];
useFMRI         = true;

% TODO: if the button was not pressed is missed the correct answer or should it be excluded from the data?

for trial = 1:length(log.data.Condition)
    onset(end+1)        = log.data.VideoStart(trial)-log.data.VideoStart(1);
    trialType{end+1}    = char(log.data.Condition(trial));
    RatingStartTime     = log.data.Rating_vbl (trial) - 1;  % Some slack time in case subs pressed a button befor the video ended
    RatingEndTime       = log.data.Rating_times (trial) + 1; % Some slack time in case subs pressed a button after response time ended
    ratingOnset(end+1)  = log.data.Rating_vbl (trial)-log.data.VideoStart(1);
    if trial < length(log.data.Condition)
        duration(end+1)         = log.data.VideoStart(trial+1)-log.data.VideoStart(trial);
        ratingDuration(end+1)   = log.data.VideoStart(trial+1)-log.data.Rating_times(trial);
    else
        duration(end+1)         = 16;
        ratingDuration(end+1)   = 2;
    end
    
    keyPressRating      = find(log.data.timeUp >= (RatingStartTime) & log.data.timeUp <= RatingEndTime);
    % PRG: This would be a miss right?
    % make sure to first extract only the valid key presses.
    
    if useFMRI
        tmpkey              = log.data.idDown(keyPressRating); % get numbers
        tmpNotTriggers      = find(tmpkey ~= ptb.Keys.trg); % get indices which are not triggers; please re-code this to make it flexibel. DONE %JB
        keyPressRating      = keyPressRating(tmpNotTriggers);
    end
    % Ideally keyPressBetween should be one (we expect only ONE key press)
    % In case we have more key presses (ignoring the triggers keys) we take
    % the last one
    if ~isempty (keyPressRating)
        keyPressRating      = keyPressRating(end);
    else
        keyPressRating      = 0;
    end
    
    %JB VERSION (dependent on log.leftId/log.rightId which is defined in ptb(.key.left/right) for each kind of setUp)

    if keyPressRating
        switch log.data.idDown(keyPressRating)
            case ptb.Keys.one 
                value{end+1} = '1';
            case ptb.Keys.two 
                value{end+1} = '2';
            case ptb.Keys.three
                value{end+1} = '3';
            case ptb.Keys.four
                value{end+1} = '4';
            case ptb.Keys.five
                value{end+1} = '5';
            otherwise
                if fillMisses == true
                    value{end+1} = '4';
                else
                    value{end+1} = 'n/a';
                end
        end
    else
        if fillMisses == true
            value{end+1} = '4';
        else
            value{end+1} = 'n/a';
        end
    end
    % pointer/cue := number indicating image to memorize.
    % By design this number always points to the solution.
    % But then during drawing the locationOfTest
    % matrix shuffles the order of test images.
end