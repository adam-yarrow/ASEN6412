%{
    This script runs the inverse UQ problem to estimate what the actual
    value of our Cm_alpha is, given a single noisy trajectory.

    Focuses on using Bayesian estimation to extract constant and linear
    CM_alpha models.
    
%}
Startup();

%% Configuration
% Fix rng
rng(11112026,'twister');

%% Generate Trajectory Instance
% Real Traj Options
mode_reality = 'linear-cmAlpha'; % Options = 'static-cmAlpha', 'linear-cmAlpha', 'Cmalpha-GP-staticMean', 'Cmalpha-GP-movMean'

% Running Model
[yMeas, t, nTimes, CM_alphaTrue, constReality, controller] = ...
    SetupInverseUQ(mode_reality);

%% Direct Bayes Approach - Estimating Cm_alpha Only
%{
    Considers the case of searching for Cm_alpha only, even if model is of
    a higher order
%}
CmalphaMin = 1.2*constReality.aero.CM_alpha;%-2.15; 
CmalphaMax = 0.8*constReality.aero.CM_alpha; %-1.9;

% Variable of interest
nCmalpha = 500; % 500
rawCmAlpha = linspace(CmalphaMin, CmalphaMax, nCmalpha)';

mode_estimator = 'static-cmAlpha'; % Options = 'static-cmAlpha', 'linear-cmAlpha', 'Cmalpha-GP-staticMean', 'Cmalpha-GP-movMean'
covInflationSF = 1.0; % Used for GP only
fPlot = true;
nWorkers = 8;

[fTheta_Y, stats, rawData] = CmAlphaEstimation(yMeas, rawCmAlpha, mode_estimator, mode_reality, ...
                                covInflationSF, fPlot, controller, nWorkers);


%% Direct Bayes Approach - Cm_Alpha and D/Dt Cm_alpha (linear model)
CmalphaMin = 1.2*constReality.aero.CM_alpha;
CmalphaMax = 0.8*constReality.aero.CM_alpha;

% Variables of interest
nCmalpha = 500;
rawCmAlpha = linspace(CmalphaMin, CmalphaMax, nCmalpha)';
dCmAlpha = rawCmAlpha(2) - rawCmAlpha(1);

derivCmAlphaMin = -constReality.aero.CM_alpha*0.1/10; %0; % assume you know thermal loading bends you in one direction
derivCmAlphaMax = -constReality.aero.CM_alpha*0.3/10; %-constReality.aero.CM_alpha*0.4/10; % max of +4% per second
nDeriv = 100;

rawDeriv = linspace(derivCmAlphaMin, derivCmAlphaMax, nDeriv);
dDeriv = rawDeriv(2) - rawDeriv(1);

[cmAlphaGrid, derivGrid] = meshgrid(rawCmAlpha,rawDeriv);
rawTheta = [cmAlphaGrid(:), derivGrid(:)];

thetaPaths = {'aero.CM_alpha','aero.movMeanCM_alpha.deriv'};

% Integration weights 
oneDimWeights_CmAlpha = dCmAlpha/2 * [1; 2*ones(nCmalpha-2,1); 1]; % TrapZ weights
oneDimWeights_deriv = dDeriv/2 * [1; 2*ones(nDeriv-2,1); 1];

integrationWeights = oneDimWeights_CmAlpha*oneDimWeights_deriv';
integrationWeights = integrationWeights(:);

% uniform prior
fPrior = ones(nDeriv, nCmalpha) * 1/((derivCmAlphaMax - derivCmAlphaMin)*(CmalphaMax - CmalphaMin));
fPrior = fPrior(:);

nWorkers = 10;
[fTheta_Y_2D, stats_2D, rawData_2D] = LinearCmAlphaEstimation(yMeas, ...
                                 rawTheta, 'linear-cmAlpha', mode_reality, ...
                                covInflationSF, fPlot, controller, nWorkers,...
                                integrationWeights, fPrior, thetaPaths);

%{
 TODO - could look at malhalanobis distance between true point and
 estimate.

Could also look at 2D KS test for how close the distributions match or KL
divergence?
%}
