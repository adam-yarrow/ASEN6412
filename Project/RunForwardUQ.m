%{
    Runs foward UQ to understand uncertainty in trajectory as a function of
    input parameters.

    Uses stochastic collocation to get mean and covariance of trajectories
    with respect to time.
    
    Compares the result with classical MC methods.

    Checks the convergence of the SC results.
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
Nmc = 500;

% Stochastic Collocation Params
levelSCC = 8; % Level for smolyak sparse grid


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

% Visualise trajectories
plotSimulationResults(rawData_SC.simData,[]);

%% Compare Stochastic Colloc and MC Estimates
meanAbsDiff = (xStats_SC.mean - xStats_MC.mean(:,:,end));
sigmaAbsDiff = (xStats_SC.cov .^0.5 - xStats_MC.cov(:,:,:,end) .^0.5);

alphaMeanAbsDiff = (alphaStats_SC.mean - alphaStats_MC.mean(:,end));
alphaSigmaAbsDiff = (alphaStats_SC.var .^0.5 - alphaStats_MC.var(:,end) .^0.5);

plotMCvsSC_StateMoments(const, times, meanAbsDiff, sigmaAbsDiff, ...
                               alphaMeanAbsDiff, alphaSigmaAbsDiff);  

%%  Stochastic Collocation Convergence
%{
    Treat the original stochastic collocation run as the reference and then
    see how the mean and covariance of the trajectory converges for reduced
    smolyak levels.
%}
cvgLevelsSC = 1:levelSCC-1;
comparisonTol = 1E-14;

[xStatsCvg_SC, alphaStatsCvg_SC, xStatsLogDiffCvg_SC, alphaStatsLogDiffCvg_SC] = ...
    CheckStochasticCollocConvergence(...
                                cvgLevelsSC, Y, rawData_SC.simData, ...
                                xStats_SC, alphaStats_SC,...
                                const, controller, stochasticParam, ...
                                comparisonTol);


% Plotting
cvgLevelsSCnames = cellfun(@num2str, num2cell(cvgLevelsSC), 'UniformOutput', false);

% Mean Response
plotCvgSC_Trajectories(const, times, xStatsCvg_SC, alphaStatsCvg_SC, false, ...
    cvgLevelsSCnames,'mean','mean');
plotCvgSC_Trajectories(const, times, xStatsLogDiffCvg_SC, ...
    alphaStatsLogDiffCvg_SC, true, cvgLevelsSCnames,'mean','mean');


plotCvgSC_Trajectories(const, times, xStatsCvg_SC, alphaStatsCvg_SC, false, ...
    cvgLevelsSCnames,'cov','var');
plotCvgSC_Trajectories(const, times, xStatsLogDiffCvg_SC, ...
    alphaStatsLogDiffCvg_SC, true, cvgLevelsSCnames,'cov','var');


%% TODO - finish plotting differences (vs time is probably ideal?)
stateIdx = 5;

figure();
hold on;
for i = 1:numel(cvgLevelsSC)
    plot(times, squeeze(xStatsLogDiffCvg_SC.mean(stateIdx, :,i)),'DisplayName',sprintf('q = %i',cvgLevelsSC(i)));
end
grid on;
legend();
ylabel(const.stateNames{stateIdx});


figure();
hold on;
for i = 1:numel(cvgLevelsSC)
    plot(times, squeeze(xStatsLogDiffCvg_SC.cov(stateIdx, stateIdx, :,i)),'DisplayName',sprintf('q = %i',cvgLevelsSC(i)));
end
grid on;
legend();
ylabel(sprintf('covariance diff - %s',const.stateNames{stateIdx}));


figure();
hold on;
plot(cvgLevelsSC, squeeze(xStatsLogDiffCvg_SC.mean(stateIdx, end,:)));
grid on;
legend();
ylabel(const.stateNames{stateIdx});

%% TODO - plot each mean trajectory vs time overlaid
%% TODO - plot the relative reduction in the final mean state/ covariance of the states vs time


