%% Question 1 
%{
    Simulate the rolling of two fair dice to estimate P(x+y = 7), where
    x and y are two i.i.d dies
%}
rng(1000);

maxTrials = 1E6;
n = (1:maxTrials)';

strideLength = 10;

targetSum = 7;
actualProb = 1/6;
dieMax = 6;
nDies = 2;

% Run Simulation
dieTosses = unidrnd(dieMax, maxTrials, nDies);
dieSum = sum(dieTosses, 2);
 
nMatchingTargetSum = cumsum(dieSum == targetSum);
probOfTargetSum = nMatchingTargetSum ./ n; 

% Plot
figure('Name','HW1 - Q1');
semilogx(n, probOfTargetSum);
hold on;
semilogx([1,maxTrials], [actualProb, actualProb], 'r--');
legend('Simulation Results',sprintf('True Probability'))
grid on;
xlabel('Number of Trials');
ylabel(sprintf('P(X = %d)', targetSum));
title(sprintf(...
    'Probability of Seeing A Sum of %d When Rolling Two Fair Dice vs Number of Trials',...
    targetSum));