%{
    This script runs the inverse UQ problem to estimate what the actual
    value of our Cm_alpha is, given a single noisy trajectory.

    
%}
%% Configuration
% Fix rng
rng(11112026,'twister');

% Model Baseline
Startup();
const = ModelParams(); % Baseline Parameters
controller =  @(t,x,const) AlphaController(t, x, const, @(t) getReferenceAlpha(t));

times = 0:const.dT:const.tEnd;

%% Generate Trajectory Instance
[simData] = Simulation(controller, const);


%% Generate Measurements
simDataNoisy = simData;


%% Plot Truth and Measurements
plotSimulationResults({simData, simDataNoisy},{'Truth Trajectory','Measurements'});