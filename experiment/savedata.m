function savedata(log,design,ptb)
    %.............................GET RESPONSES...........................%
    % Regardless of HOW the experiment ended.
    % Stop KbQueue data collection
    KbQueueStop(ptb.Keys.kbrd2); 
    KbQueueStop(ptb.Keys.kbrd1);     
    
    % Extract events
    if log.run ~= '3'
    while KbEventAvail(ptb.Keys.kbrd2)
        [evt, n] = KbEventGet(ptb.Keys.kbrd2);
        
%         if ptb.Keys.debug == 1
%             fprintf('Event is:\n'); disp(evt);
%             fprintf('\nNow %i events remaining.\n', n);
%         end
        
        if evt.Pressed == 1
            log.data.idDown   = [log.data.idDown; evt.Keycode];
            log.data.timeDown = [log.data.timeDown; evt.Time];
        else
            log.data.idUp   = [log.data.idUp; evt.Keycode];
            log.data.timeUp = [log.data.timeUp; evt.Time];
        end
    end
    end
    
    if strcmp(log.end,'Finished with errors')
        save([log.fileName '_design_error'],'design');
        save([log.fileName '_ptb_error'],'ptb');
        save([log.fileName '_log_error'],'log'); fprintf('\n Saved error data.... \n');
        if ~log.ETdummy; unixStr=['mv ' log.edfFile ' ' [log.edfFile] '_error.edf'];unix(unixStr);end
    elseif strcmp(log.end,'Escape')
        save([log.fileName '_design_cancelled'],'design');
        save([log.fileName '_ptb_cancelled'],'ptb');
        save([log.fileName '_log_cancelled'],'log'); fprintf('\n Saved cancelled data.... \n');
        if ~log.ETdummy; unixStr=['mv ' log.edfFile ' ' [log.edfFile] '_cancelled.edf'];unix(unixStr);end
    elseif strcmp(log.end,'Success')
        save([log.fileName '_design'],'design');
        save([log.fileName '_ptb'],'ptb');
        save([log.fileName '_log'],'log'); fprintf('\n Saved success data.... \n');
    end 
    
    if ~log.ETdummy
       
        % Shutdown Eyelink:
        Eyelink('StopRecording');
        Eyelink('CloseFile');
        WaitSecs(2);
        % get the file
%         try
%             fprintf('Receiving data file ''%s''\n',  log.edfFile);
%             status=Eyelink('ReceiveFile');
%             Waitsecs(2);
%             if status > 0
%                 fprintf('ReceiveFile status %d\n', status);
%             end
%             if 2==exist(log.edfFile, 'file')
%                 fprintf('Data file ''%s'' can be found in ''%s''\n',  log.edfFile, pwd );
%             end
%         catch rdf
%             fprintf('Problem receiving data file ''%s''\n', log.edfFile );
%             rdf;
%         end
        
        % Delete the Snd_use_oldstyle file to prevent messing up future
        % auditory experiments.
%         delete([PsychtoolboxConfigDir 'Snd_use_oldstyle.txt'])
        
    end
end
