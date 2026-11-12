%{
    This script runs the inverse UQ problem to estimate what the actual
    value of our Cm_alpha is, given a single noisy trajectory.

    
%}
%% Configuration
% Fix rng
rng(11112026,'twister');

% Model Baseline
Startup();

const = ModelParams(); 
% TODO - can put in some changes here if go for non constant etc.


controller =  @(t,x,const) AlphaController(t, x, const, @(t) getReferenceAlpha(t));

times = 0:const.dT:const.tEnd;

%% Generate Trajectory Instance
[simData] = Simulation(controller, const);
yMeas = getMeasurements(const, simData.x);

simDataNoisy = simData;
simDataNoisy.x = yMeas;

% Plot truth and measurements
plotSimulationResults({simData, simDataNoisy},{'Truth Trajectory','Measurements'});


%% Direct Bayes Approach - Estimating Static Cm_alpha

% Dummy testing fcn for now

CmalphaMin = 1.5*const.aero.CM_alpha;
CmalphaMax = 0.5*const.aero.CM_alpha;
nTheta = 100;
theta = linspace(CmalphaMin, CmalphaMax, nTheta)';

% uniform prior
fThetaPrior = ones(size(theta))*(1/(CmalphaMax - CmalphaMin));
thetaPath = {'aero.CM_alpha'};

[fTheta_Y, rawData] = DiscreteBayesFilter(theta, thetaPath, yMeas,...
                                            fThetaPrior, ...
                                            controller, const);

figure();
plot(theta, fTheta_Y);
grid on;
xlabel('Cm_alpha');


%% TODO - test this with Cmalpha being a stochastic process
%% OR - test with larger uncertainty/lack of observability?