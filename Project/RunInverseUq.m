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
CmalphaMin = 1.5*constReality.aero.CM_alpha;
CmalphaMax = 0.5*constReality.aero.CM_alpha;

% Variable of interest
nCmalpha = 500;
rawCmAlpha = linspace(CmalphaMin, CmalphaMax, nCmalpha)';
dCmAlpha = rawCmAlpha(2) - rawCmAlpha(1);

mode_estimator = 'static-cmAlpha'; % Options = 'static-cmAlpha', 'linear-cmAlpha', 'Cmalpha-GP-staticMean', 'Cmalpha-GP-movMean'
covInflationSF = 0.5; % Used for GP only
fPlot = true;
nWorkers = 8;

[fTheta_Y, stats, rawData] = CmAlphaEstimation(yMeas, rawCmAlpha, mode_estimator, mode_reality, ...
                                covInflationSF, fPlot, controller, nWorkers);


%% Direct Bayes Approach - Cm_Alpha and D/Dt Cm_alpha (linear model)
CmalphaMin = 1.1*constReality.aero.CM_alpha;
CmalphaMax = 0.9*constReality.aero.CM_alpha;

% Variables of interest
nCmalpha = 50;
rawCmAlpha = linspace(CmalphaMin, CmalphaMax, nCmalpha)';
dCmAlpha = rawCmAlpha(2) - rawCmAlpha(1);

derivCmAlphaMin = 0; % assume you know thermal loading bends you in one direction
derivCmAlphaMax = -constReality.aero.CM_alpha*0.4/10; % max of +4% per second
nDeriv = 50;

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

nWorkers = 6;
[fTheta_Y_2D, stats_2D, rawData_2D] = LinearCmAlphaEstimation(yMeas, ...
                                 rawTheta, 'linear-cmAlpha', mode_reality, ...
                                covInflationSF, fPlot, controller, nWorkers,...
                                integrationWeights, fPrior, thetaPaths);

%% TODO - could do p-test for each different case to see which is the most likely?

%% TODO - Cm_alpha and d/dt of Cmalpha grid


%% NLS Approach Test
covInflationSF = 1;
theta0 = [-1; 0; 0;  0]; % 3rd order poly
mode_estimator = 'poly-cmAlpha';

[theta_NLS, P_NLS, nlsDetails] = RunNLS(mode_estimator, theta0, covInflationSF, yMeas);

% Only works for linear model
if strcmpi(mode_reality,'linear-cmAlpha')
    % Plot Solver outputs
    thetaTrue = [const.aero.CM_alpha; const.aero.movMeanCM_alpha.deriv];
    plotNLS(nlsDetails, thetaTrue);
    
    % Error In Coeffs Est
    thetaRelErrorPercent = (theta_NLS - thetaTrue) ./ thetaTrue * 100;
end

if strcmpi(mode_reality,'Cmalpha-GP-staticMean')
    ReconstructCM_AlphaGP(constReality.aero.GP.Y, simData.times', constReality.aero.GP.dataKLE);
    hold on;
    t = simData.times';
    nTimes = numel(simData.times); 
    tMatrix = zeros(nTimes, numel(theta_NLS));
    for i = 1:numel(theta_NLS)
        tMatrix(:,i) = t.^(i-1);
    end

    plot(t, tMatrix* theta_NLS);
    legend('True GP','Polynomial NLS Fit');
end


%% NLS IC Sensitivity
theta0_ConstRange = linspace(-1.5,-2.5, 10);
nTerms = 5;
polyOrder = nTerms - 1;
nIC = numel(theta0_ConstRange);

thetaEstStorage = NaN(nTerms, nIC);
PestStorage = NaN(nTerms, nTerms, nIC);
JnlsStorage = NaN(nIC);

mode_estimator = 'poly-cmAlpha';
WaitMessage = parfor_wait(nIC, 'Waitbar', true);
parfor (iThetaIC = 1:nIC, nWorkers)
    theta0 = [theta0_ConstRange(iThetaIC); zeros(nTerms-1,1)];
    [thetaEstStorage(:,iThetaIC), PestStorage(:,:,iThetaIC), nlsDetails] = ...
                    RunNLS(mode_estimator, theta0, ...
                        covInflationSF, yMeas, controller);

    JnlsStorage(iThetaIC) = nlsDetails.Jnls(end);

    WaitMessage.Send;
end

WaitMessage.Destroy

figure();
plot(theta0_ConstRange, JnlsStorage);
grid on;
xlabel('CM_alpha Constant Param IC');
ylabel('Jnls');
title('Sensitivity of Jnls to IC');

% Plot diff realisations
t = simData.times';
tMatrix = BuildCmAlphaSensitivityMatrix(t, nTerms);

CM_alphaTrue = ReconstructCM_AlphaGP(constReality.aero.GP.Y, t, constReality.aero.GP.dataKLE);
hold on;
displayNames = cell(nIC + 1, 1);
displayNames{1} = 'True GP';
MSE_NLS_IC = NaN(nIC,1);
for i = 1:nIC
    CmAlpha_NLS = tMatrix * thetaEstStorage(:,i);
    MSE_NLS_IC(i) = sum((CmAlpha_NLS - CM_alphaTrue).^2);
    plot(t, CmAlpha_NLS,'--');
    displayNames{i+1} = sprintf('IC = %3.3f',theta0_ConstRange(i));
end
legend(displayNames);

% MSE vs IC
figure();
plot(theta0_ConstRange', MSE_NLS_IC);
grid on;
xlabel('CM_alpha Constant Param IC');
ylabel('MSE in CM_{\alpha}');
title('NLS MSE vs NLS IC');

%-1.94 with J = 0.452
%% Bisection approach
maxIts = 10;
theta0RelTol = 1E-2;
[theta, P, nlsDetails] = GoldenSearchNLS(-1.5, -2.5, maxIts, theta0RelTol, ...
                    nTerms, mode_estimator, covInflationSF, yMeas, ...
                    controller, true);