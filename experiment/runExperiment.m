function runExperiment()

subNr = input('Enter subject Nr: ','s');
for b = 1:3
    for r = 1:5
        block = num2str (b);
        run = num2str (r);
        MagicExperiment (subNr, block, run)
    end
end
end