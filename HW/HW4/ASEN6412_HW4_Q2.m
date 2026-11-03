%{
    Running this script executes the code necessary to generate all plots
    for HW4 Q2.

    The library path includes all codes from CANVAS plus my KLE solver from
    HW2.
%}
%% Startup
addpath(genpath('..\..\Library\'));

% Fix random seed for repeatability
rng(01112026);

%% Fixed Parameters
% G(x,w) 
muK = 1.0;
sigmaKxx = 0.1;
lKxx = 2;

% Truncation Terms
d = 2;

% Domain
xMax = 1;
xMin = 0;
Nx = 101; % From Previous HW3 solution
x = linspace(xMin, xMax, Nx)'; % Fixed domain
xHalfIdx = find(x == 0.5); if isempty(xHalfIdx); error('Invalid halfway point'); end

% KLE Random Var Uniform Domain
Ya = [-1, -1];
Yb = [1, 1];

% Boundary Conditions
uBC = [0, 0]; % Corresponds to x = [0, 1]

% Stochastic Collocation - PDF scale factor
uniformScaleFactor = 1/2 * 1/2; % 2D uniform distribution scale factor


%% Setup KLE Expansion Common Parameters
% Get d eigenfunctions
lambdaKxx = NaN(d,1);
omegaKxx = NaN(d,1);
aKxx = (xMax - xMin)/2; % Calculate domain boundary

for i = 1:d
    omegaKxx(i) = getOmegaKLE(i,aKxx,lKxx);
    lambdaKxx(i) = getLambdaKLE(lKxx,omegaKxx(i),sigmaKxx);
end

% Generate KLE Data
meanKFunc = @(x) muK * ones(size(x));
dataKLE = generateKLE_Struct(lambdaKxx, omegaKxx, meanKFunc, aKxx, d);

% Generate Phi Matrix to save time for our grid
dx = x(2) - x(1);
x_itHalf = x - dx/2; 
x_itHalf = x_itHalf(2:end); % First entry outside of domain -> remove

% Pre-make eigenfunctions - I AM SPEED!
% Eigenfunctions don't change between grids
[phiMatrix] = MakePhiMatrix(x_itHalf, dataKLE);

dataKLE.x_itHalf = x_itHalf;
dataKLE.phiMatrix = phiMatrix;


%% Standard MC Solution
N_MC = [1E1, 1E2, 1E3, 1E4, 1E5];
NumMC = numel(N_MC);

uMC = cell(NumMC,1);
uStatsMC = cell(NumMC,1);

mu_MC = NaN(Nx, NumMC);
var_MC = NaN(Nx, NumMC);

for iMC = 1:NumMC
    Y_MC = Ya' + (Yb - Ya)'.*rand(2, N_MC(iMC));
    [uMC{iMC}, uStatsMC{iMC}] = MonteCarloPDE_HW4(x, uBC, N_MC(iMC), ...
                            dataKLE, Y_MC);
    mu_MC(:,iMC) = uStatsMC{iMC}.mu;
    var_MC(:,iMC) = uStatsMC{iMC}.var;
end

% Plotting
plotUxVsNpoints('MC', x, mu_MC, var_MC, N_MC, true);
plotUxVsNpointsDifference('MC', x, mu_MC(:,1:end-1), var_MC(:,1:end-1), ...
                            mu_MC(:,end), var_MC(:,end), N_MC(1:end-1), true);

%% LHS Solution - <u(x)> and var(u(x))
N_LHS = [1E1, 1E2, 1E3, 1E4, 1E5];
NumLHS = numel(N_LHS);

fPlotLHS = false;

uLHS = cell(NumLHS,1);
uStatsLHS = cell(NumLHS,1);

mu_LHS = NaN(Nx, NumLHS);
var_LHS = NaN(Nx, NumLHS);

for iLHS = 1:NumLHS
    Y_LHS = getLHS_SamplesUniformDist(N_LHS(iLHS), Ya, Yb);     
    [uLHS{iLHS}, uStatsLHS{iLHS}] = MonteCarloPDE_HW4(x, uBC, N_LHS(iLHS), ...
                            dataKLE, Y_LHS);
    mu_LHS(:,iLHS) = uStatsLHS{iLHS}.mu;
    var_LHS(:,iLHS) = uStatsLHS{iLHS}.var;

end

% Plotting
plotUxVsNpoints('LHS', x, mu_LHS, var_LHS, N_LHS, true);
plotUxVsNpointsDifference('LHS', x, mu_LHS(:,1:end-1), var_LHS(:,1:end-1), ...
                            mu_LHS(:,end), var_LHS(:,end), N_LHS(1:end-1), true);

%% Stochastic Collocation - with CC Rule
levels_CC = [0, 1, 2, 3, 4, 5];
lRefIdx_CC = 6; % Corresponds to N = 33 points
N_levels_CC = [1, 2.^levels_CC(2:end) + 1];
Nsamples_CC = N_levels_CC.^2;
NumLevels_CC = numel(levels_CC);

% Pre - Build Flattened Tensor Points
Y_Tensors_CC = cell(NumLevels_CC,1);
W_Tensor_CC = cell(NumLevels_CC,1);

for iLevel = 1:NumLevels_CC
    [yi,wi]=spquad(1,levels_CC(iLevel));

    [Y1, Y2] = meshgrid(yi, yi);

    Y_Tensors_CC{iLevel} = [Y1(:), Y2(:)];

    W_TensorTemp = wi*wi';
    W_Tensor_CC{iLevel} = W_TensorTemp(:);
end

[mu_CC, var_CC] = StochasticColloc(Y_Tensors_CC, W_Tensor_CC, Nx, uBC, ...
    dataKLE, dx, uniformScaleFactor);

% Reference Solution
muRef = mu_CC(:,lRefIdx_CC);
varRef = var_CC(:,lRefIdx_CC);

% Plotting
plotUxVsNpoints('CC - Stochastic Collocation', x, mu_CC, var_CC, Nsamples_CC, false);
plotUxVsNpointsDifference('CC - Stochastic Collocation', x, ...
                            mu_CC(:,1:end-1), var_CC(:,1:end-1), ...
                            muRef, varRef, Nsamples_CC(1:end-1), false);

%% Stochastic Collocation - Smolyak
levels_SCC = [1, 2, 3, 4, 5, 6, 7];

NumLevels_SCC = numel(levels_SCC);

Y_Tensors_SCC = cell(NumLevels_SCC,1);
W_Tensor_SCC = cell(NumLevels_SCC,1);

N_levels_SCC = NaN(NumLevels_SCC, 1);

for iLevel = 1:NumLevels_SCC
    [Y_Tensors_SCC{iLevel},W_Tensor_SCC{iLevel}] = spquad(2, ...
                                                    levels_SCC(iLevel));
    N_levels_SCC(iLevel) = numel(W_Tensor_SCC{iLevel});
end

[mu_SCC, var_SCC] = StochasticColloc(Y_Tensors_SCC, W_Tensor_SCC, Nx, uBC, ...
    dataKLE, dx, uniformScaleFactor);

% Plotting
plotUxVsNpoints('Sparse CC - Stochastic Collocation', x, mu_SCC, var_SCC, N_levels_SCC, false);
plotUxVsNpointsDifference('Sparse CC - Stochastic Collocation', x, mu_SCC, var_SCC, ...
                            muRef, varRef, N_levels_SCC, false);


%% Comparison of All Methods at x = 0.5
% Absolute Difference
muErr = getLogAbsErr(mu_MC, mu_LHS, mu_CC, mu_SCC, muRef, xHalfIdx);
varErr = getLogAbsErr(var_MC, var_LHS, var_CC, var_SCC, varRef, xHalfIdx);

% Plotting
errNames = {'MC','LHS','Clenshaw-Curtis','Smolyak Sparse Grid'};
errNsamples = {N_MC, N_LHS, Nsamples_CC, N_levels_SCC};
% errAccLevels = {[],[],levels_CC,levels_SCC};

figure('Name','Convergence Rate of Algorithms');
subplot(1,2,1);
plotErrConvergence(errNames, errNsamples, muErr, 'log_{10}(|E[u] - \mu_{ref}|)', 'Mean');
subplot(1,2,2);
plotErrConvergence(errNames, errNsamples, varErr, 'log_{10}(|Var(u) - \sigma_{ref}^2|)', 'Variance');


%% Calculate min{d} required to obtain same MSE with l = 0.1
lambdaOptions = struct();
lambdaOptions.mode = 'approximate';
lambdaOptions.k = 200; % Assuming 100 eigenvalues beyond current eval adequate
lambdaOptions.maxIts = 100;
lambdaOptions.fPlot = true;

truncErrorMSEGuess = 0.01; % Some small value to get lots of evals

[~, lambdaOrig, ~] = findMinNumLambdaKLE(truncErrorMSEGuess, ...
                                            aKxx, lKxx, sigmaKxx, lambdaOptions);


% What number of lambda do we need to achieve the same MSE approximation as
% the d = 2 case?
truncErrorExact = 1 - (lambdaOrig(1) + lambdaOrig(2)) / sum(lambdaOrig);
[dSmallL, lambdaSmallL, ~] = findMinNumLambdaKLE(truncErrorExact, ...
                                            aKxx, 0.1, sigmaKxx, lambdaOptions);

% Plot E-vals together to show impact
figure('Name','Impact of Correlation Length');
semilogy(1:numel(lambdaOrig), lambdaOrig ./ sum(lambdaOrig),'r.','MarkerSize',8);
hold on;
semilogy(1:numel(lambdaSmallL), lambdaSmallL ./ sum(lambdaSmallL),'b.','MarkerSize',8);
grid on;
xlabel('i','FontSize',14)
txt = '$\lambda_i $ / $ \sum_{j=1}^{k}\lambda_j$';
ylabel(txt,'Interpreter','latex','FontSize',14);
title('Normalized Eigenvalues vs Mode Number','FontSize',14)
legend('l = 2','l = 0.1');



%% Helper Functions
function plotErrConvergence(errNames, errNsamples, errData, yLabelTxt, titleSpecific)
    for i = 1:numel(errNames)
        plot(log10(errNsamples{i}), errData{i},'o-','DisplayName', errNames{i},...
               'LineWidth',1.5);
        hold on;    
    end
    grid on;
    xlabel('log_{10}(N) samples');
    ylabel(yLabelTxt);
    title(sprintf('Convergence of %s of u(0.5) vs Number of Samples',titleSpecific));
    legend();
end

function absErr = getLogAbsErr(mc, lhs, cc, scc, ref, xHalfIdx)
    errMC = log10(abs(mc(xHalfIdx,:) - ref(xHalfIdx))');
    errLHS = log10(abs(lhs(xHalfIdx,:) - ref(xHalfIdx))');
    errCC = log10(abs(cc(xHalfIdx,:) - ref(xHalfIdx))');
    errSCC = log10(abs(scc(xHalfIdx,:) - ref(xHalfIdx))');

    absErr = {errMC, errLHS, errCC, errSCC};
end