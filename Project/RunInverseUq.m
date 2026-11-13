%{
    This script runs the inverse UQ problem to estimate what the actual
    value of our Cm_alpha is, given a single noisy trajectory.

    
%}
%% Configuration
% Fix rng
rng(11112026,'twister');

Startup();

% Model Baseline
controller =  @(t,x,const) AlphaController(t, x, const, @(t) getReferenceAlpha(t));


%% Generate Trajectory Instance
% Real Simulation
mode_reality = 'Cmalpha-GP-staticMean'; % Options = 'static-cmAlpha', 'linear-cmAlpha', 'Cmalpha-GP-staticMean', 'Cmalpha-GP-movMean'
constReality = BuildInverseUQParams(mode_reality);

[simData] = Simulation(controller, constReality);
yMeas = getMeasurements(constReality, simData.x);

simDataNoisy = simData;
simDataNoisy.x = yMeas;

% Plot truth and measurements
plotSimulationResults({simData, simDataNoisy},{'Truth Trajectory','Measurements'});

%% Direct Bayes Approach - Estimating Cm_alpha Only
%{
    Considers the case of searching for Cm_alpha only, even if model is of
    a higher order
%}
CmalphaMin = 1.5*constEstimator.aero.CM_alpha;
CmalphaMax = 0.5*constEstimator.aero.CM_alpha;

% Variable of interest
nCmalpha = 100;
rawCmAlpha = linspace(CmalphaMin, CmalphaMax, nCmalpha)';
dCmAlpha = rawCmAlpha(2) - rawCmAlpha(1);

mode_estimator = 'Cmalpha-GP-staticMean'; % Options = 'static-cmAlpha', 'linear-cmAlpha', 'Cmalpha-GP-staticMean', 'Cmalpha-GP-movMean'
covInflationSF = 0.5; % Used for GP only
fPlot = true;
nWorkers = 6;

[fTheta_Y, stats, rawData] = CmAlphaEstimation(yMeas, rawCmAlpha, mode_estimator, mode_reality, ...
                                covInflationSF, fPlot, controller, nWorkers);


%% Direct Bayes Approach - Cm_Alpha and D/Dt Cm_alpha (linear model)
derivCmAlphaMin = 0; % assume you know thermal loading bends you in one direction
derivCmAlphaMax = constEstimator.aero.CM_alpha*0.4/10; % max of 4% per second
nDeriv = 100;

rawDeriv = linspace(derivCmAlphaMin, derivCmAlphaMax, nDeriv);
dDeriv = rawDeriv(2) - rawDeriv(1);

[derivGrid, cmAlphaGrid] = meshgrid(rawDeriv, rawCmAlpha);
rawTheta = [derivGrid(:), cmAlphaGrid(:)];


% Integration weights 
oneDimWeights_CmAlpha = dCmAlpha/2 * [1; ones(nCmAlpha-2,1); 1]; % TrapZ weights
oneDimWeights_deriv = dDeriv/2 * [1; ones(nDeriv-2,1); 1];

integrationWeights = oneDimWeights_CmAlpha*oneDimWeights';

% uniform prior
fPrior = ones(nDeriv, nCmalpha) * 1/((derivCmAlphaMax - derivCmAlphaMin)*(CmalphaMax - CmalphaMin));
fPrior = fPrior(:);

[fTheta_Y_2D, stats_2D, rawData_2D] = LinearCmAlphaEstimation(yMeas, rawTheta, mode_estimator, mode_reality, ...
                                covInflationSF, fPlot, controller, nWorkers,...
                                integrationWeights, fPrior);

%% TODO - could do p-test for each different case to see which is the most likely?

%% TODO - Cm_alpha and d/dt of Cmalpha grid