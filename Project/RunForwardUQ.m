%{
    Runs foward UQ to understand uncertainty in trajectory as a function of
    input parameters.

    Uses stochastic collocation to get mean and covariance of trajectories
    with respect to time.
    
    Compares the result with classical MC methods
%}
%% Configuration
% Fix rng
rng(11062026,'twister');

% Model Baseline
Startup();
const = ModelParams(); % Baseline Parameters
controller =  @(t,x,const) AlphaController(t, x, const, @(t) getReferenceAlpha(t));

times = 0:const.dT:const.tEnd;

% Cmq Active and static model
const.aero.dynDeriv.model = 'static'; 

% Random Inputs
[stochasticParam] = ForwardUqPerturbationVars();
dInputs = numel(stochasticParam);

% MC Params
Nmc = 100;

% Stochastic Collocation Params
levelSCC = 3; % Level for smolyak sparse grid


%% Baseline MC Method
% Run MC
[xStats_MC, alphaStats_MC, rawData_MC] = MonteCarloFwdUQ(const,Nmc,controller,...
    stochasticParam);

% Plotting State Stats vs Time
plotStateMoments(const, times, xStats_MC, alphaStats_MC, 'MC');

% Plotting State Convergence At End time
tToPlotConvergenceIdx = find(times == const.tEnd);
plotConvergenceOfStateMoments(const, tToPlotConvergenceIdx, xStats_MC,...
                                alphaStats_MC, Nmc);

% runNames = cellfun(@num2str, num2cell(1:Nmc), 'UniformOutput', false);
plotSimulationResults(rawData_MC.simData,[]);

%% Stochastic Colloc
% Want to use a sparse CC grid for implementation ease & ability to check
% mesh convergence easier.
% Generate Raw Weights and Yi values (Yi = [-1,1])
[Y, W] = spquad(dInputs, levelSCC);

% Run Fwd UQ
[xStats_SC, alphaStats_SC, rawData_SC] = StochasticColloc(const, Y, W, ...
                         controller, stochasticParam);

% Plotting State Stats vs Time
plotStateMoments(const, times, xStats_SC, alphaStats_SC, 'SC');

%% Compare Stochastic Colloc and MC Estimates

%% TODO - validate working correctly by comparing with MC results
%% TODO - most worried about the multi-dimensional weightings

% Testing
% [stochasticParam] = ForwardUqPerturbationVars();
% [Y, W] = spquad(dInputs, 3);
% [stochasticValues, rawParams, Wcorrected] = MapSCCtoDistributions(...
%                     stochasticParam, Y, W);



