function MagicGenerateCounterBalance()

number_of_participants = 36;
number_of_blocks       = 3;

design.possible_perms  = perms(1:number_of_blocks);
design.subjSequence    = repmat(1:length(design.possible_perms),1,number_of_participants);

% design.matrix.objectLevels = ["Ball_" "Ball_" "Ball_"];
design.matrix.objectLevels = ["Ball_" "Card_" "Stick_"];
design.matrix.trickLevels  = ["Appear1_" "Appear2_" "Vanish1_" "Vanish2_" "Change1_" "Change2_"];
design.matrix.effectLevels = ["Magic" "Magic" "Control"]; % We show the magic videos twice as much
design.matrix.delimiter    = '_';

[design.matrix.Object, design.matrix.Trick, design.matrix.Effect] = BalanceFactors(1,0, design.matrix.objectLevels, design.matrix.trickLevels, design.matrix.effectLevels);
design.Condition = strcat(design.matrix.Object, design.matrix.Trick, design.matrix.Effect);
design.Condition = reshape(design.Condition,[length(design.Condition)/length(design.matrix.objectLevels) length(design.matrix.objectLevels)]);
% design.Condition = [design.Condition; ["Ball_Surprise1" "Ball_Surprise1" ...
%                     "Ball_Surprise1"; "Ball_Surprise2" "Ball_Surprise2" "Ball_Surprise2"]];
design.Surprises = ["Ball_Surprise1" "Card_Surprise1" "Stick_Surprise1"; ...
    "Ball_Surprise2" "Card_Surprise2" "Stick_Surprise2";...
    "Ball_Surprise3" "Card_Surprise3" "Stick_Surprise3"];
design.Surprises = repmat (design.Surprises, [2 1]);
design.Condition = [design.Condition; design.Surprises];

[design.RevealingMatrix.Object, design.RevealingMatrix.Trick] = BalanceFactors(1,0, design.matrix.objectLevels, design.matrix.trickLevels);
design.Revealing = strcat(design.RevealingMatrix.Object, design.RevealingMatrix.Trick, "Revealing");
design.Revealing = reshape(design.Revealing,[length(design.Revealing)/length(design.matrix.objectLevels)...
    length(design.matrix.objectLevels)]);

design.Tricks = strcat(design.RevealingMatrix.Object, design.RevealingMatrix.Trick, "Magic");
design.Tricks = reshape(design.Tricks,[length(design.Tricks)/length(design.matrix.objectLevels)...
    length(design.matrix.objectLevels)]);

design.Tricks = repmat(design.Tricks,1,1,number_of_participants);

design.Revealing = repmat(design.Revealing,1,1,number_of_participants);
design.Condition = repmat(design.Condition,1,1,number_of_participants);

for i = 1:number_of_participants
    design.Tricks(:,:,i)    = design.Tricks(:,design.possible_perms(design.subjSequence(i),:),i);
    design.Condition(:,:,i) = design.Condition(:,design.possible_perms(design.subjSequence(i),:),i);
    design.Revealing(:,:,i) = design.Revealing(:,design.possible_perms(design.subjSequence(i),:),i);
end

save('CounterBalancedSubjects','design')

end

