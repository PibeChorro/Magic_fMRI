function [ptb, log] = MagicGetTrigger(ptb,log,dummies)
%............................ TRIGGER ....................................%
% Wait for scanner triggers before starting the experiment.
% PRG 12/2015

waiting4dummies          = 1;
log.exp.dummies.nScanned = 0;
scannerTrg               = 1;

fprintf('\n========================================');
fprintf('\nWaiting for dummy triggers!');
fprintf('\n========================================');

% Wait for trigger and dummies
while waiting4dummies
    
    %if log.Keys.break; break; end
    
    %     if ptb.str.instruction % show start screen while waiting for triggers
    %         Screen('FillRect', ptb.w.id, ptb.w.bg);
    %         DrawFormattedText(ptb.w.id,[ptb.str.instruction '\n \n' ptb.str.trgstr], 'center', 'center', [0.85 0.85 0]); % GELB
    %         Screen('Flip', ptb.w.id);
    %     end
    
    if ptb.usbTrg
        %KbQueueCreate(); % for subject's responses
        KbQueueStart(ptb.Keys.kbrd2);
        %[~, firstPress] = KbQueueCheck(pys.kbrd2);
        
        [~, firstPress] = KbQueueCheck(ptb.Keys.kbrd2);
        
        % If trigger arrived
        if firstPress(ptb.Keys.trg)
            
            log.Keys.trig(scannerTrg) = firstPress(ptb.Keys.trg);
            
            
            if scannerTrg == 1
                log.exp.dummies.tFirstDummyArrived = firstPress(ptb.Keys.trg);
            end
            
            if log.exp.dummies.nScanned < dummies
                log.exp.dummies.nScanned = log.exp.dummies.nScanned+1;
                fprintf('\nGot dummy %d of %d', log.exp.dummies.nScanned, dummies);
                log.exp.dummies.t(scannerTrg) = firstPress(ptb.Keys.trg);
                scannerTrg = scannerTrg + 1;
            else
                
                fprintf('\n');
                fprintf('\n========================================');
                fprintf('\nStarting stimulus ...');
                fprintf('\n========================================\n');
                
                ptb.Keys.FirstVolOfExp = firstPress(ptb.Keys.trg);
                waiting4dummies       = 0;
                
                log.exp.dummies.t     = log.exp.dummies.t';
                log.Keys.trig         = log.Keys.trig';
            end
        end
    end
    

end
