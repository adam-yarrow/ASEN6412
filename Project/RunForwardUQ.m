%{
    Runs foward UQ to understand uncertainty in trajectory as a function of
    input parameters
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

% MC Params
Nmc = 25;
[pathsToPerturb, stochasticParam] = ForwardUqPerturbationVars();

% Stochastic Collocation Params


%% Baseline MC Method
% Run MC
[xStats_MC, alphaStats_MC, rawData_MC] = MonteCarloFwdUQ(const,Nmc,controller,...
    pathsToPerturb,stochasticParam);

% Plotting State Stats vs Time
plotStateMoments(const, times, xStats_MC, alphaStats_MC);

% Plotting State Convergence At End time
tToPlotConvergenceIdx = find(times == const.tEnd);
plotConvergenceOfStateMoments(const, tToPlotConvergenceIdx, xStats_MC,...
                                alphaStats_MC, Nmc);

% runNames = cellfun(@num2str, num2cell(1:Nmc), 'UniformOutput', false);
plotSimulationResults(rawData_MC.simData,[]);

%TODO - plot convergence of final state for variables of interest and
%variance - can i plot 3D confidence ellipses over the trajectory?

%% Stochastic Colloc

