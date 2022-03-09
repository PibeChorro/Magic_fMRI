function cleanup(log)

    % Restore keyboard output to Matlab:
    ListenChar(0);

    if ~log.ETdummy
        Eyelink('Shutdown');
    end
    
    % Close window:
    sca;
end