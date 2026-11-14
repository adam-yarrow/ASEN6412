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
nTimes = size(yMeas,2);
t = simData.times';

simDataNoisy = simData;
simDataNoisy.x = yMeas;

% Plot truth and measurements
plotSimulationResults({simData, simDataNoisy},{'Truth Trajectory','Measurements'});

CM_alphaTrue = ReconstructCM_AlphaGP(constReality.aero.GP.Y, t, ...
                constReality.aero.GP.dataKLE);


%% Polynomial Config
maxItsGoldenSearch = 10;
CmAlphaMin = -2.5;
CmAlphaMax = -1.5;

theta0RelTolGoldenSearch = 1E-2;

mode_estimator = 'poly-cmAlpha';

covInflationSF = 1;

%% Find Best Polynomial Order Model
polyOrder = [3,4,5];
nTestPoly = numel(polyOrder);

thetaVsPoly = cell(nTestPoly,1);
PVsPoly = cell(nTestPoly,1);
JvsPoly = NaN(nTestPoly,1);

CmAlphaVsPoly = NaN(nTimes, nTestPoly);
MSE_VsPoly = NaN(nTestPoly,1);

% WaitMessage = parfor_wait(nTestPoly, 'Waitbar', true);
runTimes = NaN(nTestPoly,1);
for iOrder = 1:nTestPoly
    tic;
    nTerms = polyOrder(iOrder) + 1;
    [thetaVsPoly{iOrder}, PVsPoly{iOrder}, nlsDetails] = ...
                GoldenSearchNLS(CmAlphaMin, CmAlphaMax, maxItsGoldenSearch,...
                    theta0RelTolGoldenSearch, ...
                    nTerms, mode_estimator, covInflationSF, yMeas, ...
                    controller, false);

    JvsPoly(iOrder) = nlsDetails.Jnls(end);

    tMatrix = BuildCmAlphaSensitivityMatrix(simData.times', nTerms);
    CmAlphaVsPoly(:,iOrder) = tMatrix * thetaVsPoly{iOrder};
    MSE_VsPoly(iOrder) = sum((CmAlphaVsPoly(:,iOrder) - CM_alphaTrue).^2);

    runTimes(iOrder) = toc;

    % WaitMessage.Send;
end

% WaitMessage.Destroy;


% Plotting
figure();
plot(polyOrder, JvsPoly);
grid on;

ReconstructCM_AlphaGP(constReality.aero.GP.Y, t, ...
                constReality.aero.GP.dataKLE)
dispNames = cell(nTestPoly+1,1);
dispNames{1} = 'True CM_{\alpha| - GP';
for i = 1:nTestPoly
    plot(t, CmAlphaVsPoly(:,i),'--');
    dispNames{i+1} = sprintf('Polynomial Approximant Order = %i', polyOrder(i));
end
legend(dispNames);
