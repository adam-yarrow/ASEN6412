%% Checking orthonormality of Polynomials
maxOrder = 3;

figure();
hold on;
y = -1:0.1:1;
for i = 0:maxOrder
    plot(y, LegendrePoly(y, i),'DisplayName',num2str(i));        
end
grid on;
xlabel('y');
legend();

%% Assuming you ran Fwd UQ and got MC results first

%% Basis reconstruction
times = 0:const.dT:const.tEnd; times = times';
nTimes = numel(times);
d = dInputs; % TODO - will need to pull this from data first and include Y if GP
N = numel(rawData_MC.simData); % Number of simulations available
p = 5; % Total order basis
P = factorial(p+d)/(factorial(p)*factorial(d)) - 1;
psiMatrix = NaN(N, P+1);
%C = NaN(P+1, nTimes);
% Solve psiMatrix*C = uMatrix to get ci(t) coefficients for each basis
% function

%% Convert MC data into useful formats
% Get Input Vectors for each simulation - assumes they are numerical values
y = [rawData_MC.rawStochasticValues{:}]'; % d x N (N simulations)

% Assuming uniform distribution - Hence Legendre 
%% NOTE: if use uniform basis for GP then this all needs to change...
yNorm = NaN(size(y));
for iDim = 1:d
    b = stochasticParam{iDim}.max;
    a = stochasticParam{iDim}.min;

    yNorm(iDim,:) =  2*(y(iDim,:) - (b+a)/2)/(b-a);
end

% Re-order raw states to be: N x nTimes x nStates
uMatrix = permute(rawData_MC.x, [3,2,1]);

%% Build psiMatrix once for all states
[validBasisCombos] = GetValidPC_BasisOptions(d,p); % Need this for picking out correct bases for SOBOL indices
for iMC = 1:N
    psiMatrix(iMC,:) = BuildPC_LegendreBasis(yNorm(:,iMC), validBasisCombos);
end

%% Solve System for each state
Cmatrix = NaN(P+1, nTimes, const.nStates);
RMSE_PCE = NaN(const.nStates, N);
xEst = NaN(N, nTimes, const.nStates);
for iState = 1:const.nStates
    Cmatrix(:,:,iState) = psiMatrix \ uMatrix(:,:,iState);
    
    % MSE calculation
    xEst(:,:,iState) =  psiMatrix*Cmatrix(:,:,iState);
    RMSE_PCE(iState,:) = sqrt(sum((xEst(:,:,iState) - uMatrix(:,:,iState)).^2, 2));
end

%% Plotting - Checking how good the solution is
figure('Name','PCE RMSE Histograms');
for iState = 1:const.nStates
    subplot(2,3,iState);
    histogram(RMSE_PCE(iState,:)*const.statePlottingSF(iState));
    grid on;
    xlabel(sprintf('%s (%s)',const.stateNames{iState}, const.statePlottingUnits{iState}))
    ylabel('Counts');
end
sgtitle('PCE RMSE Histograms');

%% TODO - generate fits for all states and look at MSE over all samples 
% e.g. plot MSE in state vs run number and compare if acceptable
