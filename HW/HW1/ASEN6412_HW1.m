%% Question 1 
%{
    Simulate the rolling of two fair dice to estimate P(x+y = 7), where
    x and y are two i.i.d dies
%}
%% TODO - replace n with a logspace and go out longer
n = [50,100, 1000, 10000];

nRuns = numel(n);

targetSum = 7;
dieMax = 6;
nDies = 2;

probOfTargetSum = zeros(nRuns,1);

% Run Simulation
runCounter = 1;

%% TODO - change this loop to just simulate once and then truncate the
% results.... to plot different n

for nCurrent = n
    % Toss the Die
    dieTosses = unidrnd(dieMax, nCurrent, nDies);
    dieSum = sum(dieTosses, 2);
    
    % Number matching target
    nMatchingTargetSum = sum(dieSum == targetSum);

    % Results
    probOfTargetSum(runCounter) = nMatchingTargetSum / nCurrent;

    runCounter = runCounter + 1;    
end

% Plot
figure('Name','HW1 - Q1');
hold on;
plot(n, probOfTargetSum);
yline(1/6, 'r--');
legend('Simulation Results','True Probability')
grid on;
xlabel('Number of Trials');
ylabel(sprintf('P(x = %d', targetSum));
title('Probability of Seeing A Sum of %d When Rolling Two Fair Dice vs Trials');