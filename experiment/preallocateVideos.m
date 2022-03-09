function [design,MovieFiles] =preallocateVideos(design)

tmpMovieFile = strcat([pwd,'/../Stimuli/'], design.Condition(:),'.mp4');
for i = 1:length (tmpMovieFile)
    MovieFiles{i} = char(tmpMovieFile(i,:));
%     tmpMoviePre(i) = Screen ('OpenMovie', ptb.window, MovieFiles{i});
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


