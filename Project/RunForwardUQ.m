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

% Cmq Active and static model
const.aero.dynDeriv.model = 'static'; 

% MC Params
Nmc = 100;
[pathsToPerturb, stochasticParam] = ForwardUqPerturbationVars();

% Stochastic Collocation Params


%% Baseline MC Method
% Storage
nTimes = const.tEnd/const.dT + 1;
xStorage = NaN(const.nStates, nTimes, Nmc);
uStorage = NaN(nTimes, Nmc);
alphaStorage = NaN(nTimes, Nmc);
times = 0:const.dT:const.tEnd;

% Generate Random Inputs for each MC case
[stochasticValues, rawStochasticValues] = DrawRandomVariables(stochasticParam, Nmc);

t0 = tic;
wb = progressBar(1, Nmc, t0, [], 'MC Runs');
simDataStorage = cell(Nmc, 1);
for iMC = 1:Nmc
    % Overwrite perturbations
    constTemp = updateStruct(const, pathsToPerturb, stochasticValues{iMC});

    % Run Model
    simData = Simulation(controller, constTemp);
    simDataStorage{iMC} = simData;
    xStorage(:,:,iMC) = simData.x;
    uStorage(:,iMC) = simData.u;
    alphaStorage(:,iMC) = calcAlpha(simData.x);

    % Progress bar
    progressBar(iMC, Nmc, t0, wb);
end

% Calculate Statistics
[xStats, alphaStats] = CalculateMC_Statistics(xStorage, alphaStorage);

% Plotting State Stats vs Time
plotStateMoments(const, times, xStats, alphaStats);

% Plotting State Convergence At End time
tToPlotIdx = find(times == const.tEnd);
plotConvergenceOfStateMoments(const, tToPlotIdx, xStats, alphaStats, Nmc);

%TODO - plot convergence of final state for variables of interest and
%variance - can i plot 3D confidence ellipses over the trajectory?



%% Stochastic Colloc

