%{
    This script searches for the best polynomial approximation to the
    CM_alpha of a given trajectory as a function of time.

    It is ideally ran with the real simulation using a gaussian process, to
    help replicate model miss-match and how to handle the uncertainty
    associated with that.

    To find the best polynomial fit this script does the following:
        1. For a fixed polynomial order, use golden search to find the
        initial conditions for the zero'th order polynomial with the best
        cost function.
        2. Repeat for different polynomials.

    The rationale for searching on the zeroth order polynomial coefficient
    as the NLS initial guess as it is assumed that CM_alpha is constant,
    however, in reality it may deviate from this constant in a "slow"
    manner, hence getting a good initial guess for CM_alpha_0*t^0 should 
    help to best place the solution.

    TODO: Alternative solution would be to random search the entire
    coefficient space with MC methods, but this would be expensive and no
    better guarantees than current approach.
    
%}

Startup();

%% General Configuration
% Fix rng
rng(11112026,'twister');

%% Generate Real Trajectory Instance
% Real Traj Options
mode_reality = 'Cmalpha-GP-staticMean'; % Options = 'static-cmAlpha', 'linear-cmAlpha', 'Cmalpha-GP-staticMean', 'Cmalpha-GP-movMean'

% Running Model
[yMeas, t, nTimes, CM_alphaTrue, constReality, controller] = ...
    SetupInverseUQ(mode_reality);

%% NLS Polynomial Config
% Golden Search Opts for best initial guess for NLS
maxItsGoldenSearch = 10; % 
CmAlphaMin = -2.5; % lower bound on search
CmAlphaMax = -1.5; % upper bound on search
theta0RelTolGoldenSearch = 1E-2; % 1% tolerance on CM_alpha_0

% "Fixed" Constants
mode_estimator = 'poly-cmAlpha'; % Essentially fixed
covInflationSF = 1; %NOTE: Unused

% Poly order
polyOrder = [3,4,5]; % checking order 0 and order 1 is easier (and more 
% reliable) with discrete bayes filter, so skipping

%% Find Best Polynomial Order Model
nTestPoly = numel(polyOrder);

thetaVsPoly = cell(nTestPoly,1);
PVsPoly = cell(nTestPoly,1);
JvsPoly = NaN(nTestPoly,1);

CmAlphaVsPoly = NaN(nTimes, nTestPoly);
MSE_VsPoly = NaN(nTestPoly,1);
runTimes = NaN(nTestPoly,1);
nlsDetails = cell(nTestPoly,1);

% WaitMessage = parfor_wait(nTestPoly, 'Waitbar', true); % I was using
% par-for here but it broke a bit with a nested parfor inside goldensearch
for iOrder = 1:nTestPoly
    tic;
    nTerms = polyOrder(iOrder) + 1;
    [thetaVsPoly{iOrder}, PVsPoly{iOrder}, nlsDetails{iOrder}] = ...
                GoldenSearchNLS(CmAlphaMin, CmAlphaMax, maxItsGoldenSearch,...
                    theta0RelTolGoldenSearch, ...
                    nTerms, mode_estimator, covInflationSF, yMeas, ...
                    controller, true);

    JvsPoly(iOrder) = nlsDetails{iOrder}.Jnls(end);

    tMatrix = BuildCmAlphaSensitivityMatrix(t, nTerms);
    CmAlphaVsPoly(:,iOrder) = tMatrix * thetaVsPoly{iOrder};
    MSE_VsPoly(iOrder) = sum((CmAlphaVsPoly(:,iOrder) - CM_alphaTrue).^2);

    runTimes(iOrder) = toc;

    % WaitMessage.Send;
end
% WaitMessage.Destroy;


%% Plotting
% Cost Function
figure('Name','Inverse UQ - NLS Cost Function vs Polynomial Order');
plot(polyOrder, JvsPoly);
xlabel('Polynomial Order');
ylabel('J_{NLS}');
title('NLS Cost Function Value Vs Polynomial Fit Order for Cm_{\alpha}(t)')
grid on;

% MSE
figure('Name','Inverse UQ - MSE vs Polynomial Order');
plot(polyOrder, MSE_VsPoly);
xlabel('Polynomial Order');
ylabel('MSE');
title('MSE To True Cm_{\alpha}(t) Vs Polynomial Fit Order for Cm_{\alpha}(t)')
grid on;

% Plotting Different Realizations of CM_alpha(t)
figure('Name','Inverse UQ - CM_alpha vs Time for different Polys')
hold on;
plot(t, CM_alphaTrue,'k');
dispNames = cell(nTestPoly+1,1);
dispNames{1} = 'True C_M_{\alpha} - GP';
for i = 1:nTestPoly
    plot(t, CmAlphaVsPoly(:,i),'--');
    dispNames{i+1} = sprintf('NLS - Polynomial Order = %i', polyOrder(i));
end
legend(dispNames);
grid on;
xlabel('Time (s)');
ylabel('C_M_{\alpha}');
title('Inverse UQ - C_M_{\alpha} vs Time');
