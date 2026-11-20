% 6412 - HW5 - Q2
% Austin Bodin
% Adam Yarrow

%% Problem 2

addpath(genpath(pwd));
addpath(genpath('..\..\Library'));

load("u_samples.mat");
load("y_samples.mat");

p = 3;
d = 20;
%% Run L1 regularisation
nRandCases = 10;

% Non-Random Sampling
[meanU_Orig, varU_Orig, Ncases] = runL1Sampling(false, y_samples, u_samples, p, d);

% Run a few random samples
randMeanU = NaN(numel(Ncases), nRandCases);
randVarU = NaN(numel(Ncases), nRandCases);

for i = 1:nRandCases
    [randMeanU(:,i), randVarU(:,i), ~] = runL1Sampling(true, y_samples, u_samples, p, d);
end

%% Sample statistics
muTrue = mean(u_samples);
varTrue = var(u_samples);

%% Plot
fPlotDecrement = true;
if fPlotDecrement
    dataMean = log10(abs((meanU_Orig - muTrue) / muTrue));
    dataMeanRandom = log10(abs((randMeanU - muTrue) / muTrue));
    
    dataVar = log10(abs((varU_Orig - varTrue) / varTrue));
    dataVarRandom = log10(abs((randVarU - varTrue) / varTrue));

    yMeanName = 'E[u(y)] Log10 Relative Decrement';
    yVarName = 'Var(u(y)) Log10 Relative Decrement';
    titleModifier = 'Log Decrement ';
else
    dataMean = meanU_Orig;
    dataVar = varU_Orig;

    dataMeanRandom = randMeanU;
    dataVarRandom = randVarU;

    yMeanName = 'E[u(y)]';
    yVarName = 'Var(u(y))';
    titleModifier = '';
end

figure('Name','ASEN6412 - HW5 - Q2');
subplot(1,2,1);
hold on;
plot(Ncases, dataMean,'k','DisplayName','Samples 1:N');
for i = 1:nRandCases
    plot(Ncases, dataMeanRandom(:,i),'--','DisplayName',sprintf('Random samples Case %i', i));    
end
grid on;
legend();
xlabel('Number of Samples');
ylabel(yMeanName);
title(sprintf('%sE[u(y)] vs MC Samples - l_1 Minimization',titleModifier))

subplot(1,2,2);
hold on;
plot(Ncases, dataVar,'k','DisplayName','Samples 1:N');
for i = 1:nRandCases
    plot(Ncases, dataVarRandom(:,i),'--','DisplayName',sprintf('Random samples Case %i', i));    
end
grid on;
xlabel('Number of Samples');
ylabel(yVarName);
title(sprintf('%sVar(u(y)) vs Samples - l_1 Minimization',titleModifier))


%% Plot sparsity of problem
P = factorial(p+d) / (factorial(p) * factorial(d)) - 1;
N = numel(u_samples); %(P+1);
psi = zeros(N, P+1);
index_pc = nD_polynomial_array(d,p);
for i = 1:N        
    psi(i,:) = piset(y_samples(i,:), index_pc);
end

c = pinv(psi)*u_samples(1:end);

figure('Name','Coefficients Plot');
loglog(1:(P+1), abs(c),'o');
grid on;
xlabel('Coefficient Index');
ylabel('|c_i|');
title('Least Squares Coefficients vs Index Number')



%% Functions

function [meanU, varU, Ncases] = runL1Sampling(fRandomSampling, y_samples, u_samples, p, d)

P = factorial(p+d) / (factorial(p) * factorial(d)) - 1;
Ncases = ceil([0.2; 0.5; 0.9] * (P+1));
numN = numel(Ncases);

index_pc = nD_polynomial_array(d,p);
c_ell1 = cell(numN,1);

meanU = NaN(numN,1);
varU = NaN(numN,1);

opts = spgSetParms('iterations',6000,'verbosity',1,...
     'optTol',1e-9,'bpTol',1e-9);

sigma = 0.005;

for iN = 1:numN
    N = Ncases(iN);
    
    if fRandomSampling
        idx = randsample(1:numel(u_samples),N,false);
        y = y_samples(idx,:);
        u = u_samples(idx);
    else
        y = y_samples(1:N,:);
        u = u_samples(1:N);
    end
    % Make psi
    psi = zeros(N, P+1);
    for i = 1:N        
        psi(i,:) = piset(y(i,:), index_pc);
    end

    delta = sigma*norm(u);

    % solve system   
    [c_ell1{iN},res,gap,info] = spg_bpdn(psi,u,delta,opts);

    % extract Mean and variance
    meanU(iN) = c_ell1{iN}(1); % First order coefficient

    % note: <psi^2(y)> = 1 for legendre polynomials that have been
    % orthonomalized
    varU(iN) = sum(c_ell1{iN}(2:end).^2); % Remove first term in summation to remove mean^2
    
end
end