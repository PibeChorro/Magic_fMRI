function ShowVideo(ptb, tmpMovieFile, tmpCondition,run)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ischar (run)
    run = str2double(run);
end

tmpCondition = char(tmpCondition);
NeedForFlip = false;
% Half of the videos shall be presented mirrored along the horizontal axis.
% half of the magic videos have an 'F' for flip
% the control Videos are presented only once per run, so they can be
% flipped in the even runs
% --> so check whether the Condition has an 'F' AND we ar in an odd run OR 
% (exclusive or) we are in an even run AND it is NOT a flip condition
% in odd runs only the 'F' condition videos are mirrored and in even trial
% % the non 'F' condition videos and the control videos are flipped
if xor (mod(run,2)&&(tmpCondition(end)=='F'), (~mod(run,2)&& tmpCondition(end)~='F'))
    NeedForFlip = true;
    Screen ('glTranslate', ptb.window, ptb.xCenter, ptb.yCenter, 0);
    Screen('glScale', ptb.window, -1, 1);
    Screen ('glTranslate', ptb.window, -ptb.xCenter, -ptb.yCenter , 0);
end

% If it is a flip condition the 'F' and the '_' must be removed
if tmpCondition (end) == 'F'
    tmpCondition (end) = '';
    tmpCondition (end) = '';
end

% tmpMovieFile = strcat(pwd, '/../Stimuli/', tmpCondition, '.mp4');
% tmpMoviePre  = Screen ('OpenMovie', ptb.window, tmpMovieFile);
Screen ('PlayMovie',tmpMovieFile, 1, [], 0);

%%%%%%%%%%%%%%% Implement pimer, showing a 'Z' for magic videos%%%%%%%%%%%%
if contains (tmpCondition, 'Magic')
    Prime = 'M';
else
    Prime = 'X';
end
primeTimestampPre   = GetSecs();
PrimeOval           = [ptb.windowRect(3)/2-50,ptb.windowRect(4)/2-50 ...
    ptb.windowRect(3)/2+50,ptb.windowRect(4)/2+50];
OldFontSize = ptb.FontSize;
ptb.FontSize = Screen('TextSize', ptb.window, 60);
while true 
    % Wait for next movie frame, retrieve texture handle to it
    
    tex = Screen('GetMovieImage', ptb.window, tmpMovieFile);
    
    % Valid texture returned? A negative value means end of movie
    % reached:
    if tex<=0
        % We're done, break out of loop:
        break;
    end
    
    % Draw the new texture immediately to screen:
    Screen('DrawTexture', ptb.window, tex);
    
    % Draw Prime in the first 300 ms
    primeTimestampPost = GetSecs();
    if primeTimestampPost - primeTimestampPre < 0.5 && run~=3
        Screen ('FillOval', ptb.window, ptb.grey, PrimeOval);
        DrawFormattedText (ptb.window, Prime, 'center', 'center',ptb.black);
    end
    
    % Update display: % here make also sure you are not missing frames!
    % PRG
    Screen('Flip', ptb.window);
    % Release texture:
    Screen('Close', tex);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Not sure what exactly this part is for%
Screen('PlayMovie', tmpMovieFile, 0);   %
%Screen('CloseMovie', tmpMovieFile);     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if NeedForFlip
    Screen('glTranslate', ptb.window, ptb.xCenter, ptb.yCenter, 0);
    Screen('glScale', ptb.window, -1, 1);
    Screen('glTranslate', ptb.window, -ptb.xCenter, -ptb.yCenter , 0);
end
ptb.FontSize = Screen('TextSize', ptb.window, 30);

end

