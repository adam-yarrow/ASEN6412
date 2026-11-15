%{
    This script focuses on MC results ONLY.
        
    Objectives:
        1. Generates MC simulations of the trajectory.
        2. Generates PCE approximation of the system dynamics.
        3. Calculats the SOBOL indices of the MC solution using PCE
        approximation.    

    Assumptions:
        - All input variables are drawn from uniform distributions.
%}
Startup();

% Fix rng
rng(11142026,'twister');

%% Configuration
const = ModelParams(); % Baseline Parameters
controller =  @(t,x,const) AlphaController(t, x, const, @(t) getReferenceAlpha(t));
times = 0:const.dT:const.tEnd;

% MC Params
Nmc = 2500; 
nSamplesPerBin = 100; % For plotting interactive histograms

% PCE Parameters
maxTotalOrderPCE = 5; % (p)

% Random Variables Options: 'static-Cmalpha-wIC', 'Cmalpha-GP-staticMean'
inputVarianceMode = 'Cmalpha-GP-staticMean'; 

% Plotting
fPlotSimData = false; % Enables/disables plotting of all trajectories from MC

%% Build Random Inputs
[stochasticParam, stochasticParamGP, fGP] = ...
    ForwardUqPerturbationVars(inputVarianceMode, const);
dInputs = numel(stochasticParam);

if fGP
    [dKLE, dataKLE] = GetTruncatedKLE(stochasticParamGP);
    dInputs = dInputs + dKLE;
    stochasticParamGP.dataKLE = dataKLE;
end

%% Run MC
[xStats_MC, alphaStats_MC, rawData_MC] = MonteCarloFwdUQ(const,Nmc,controller,...
    stochasticParam, stochasticParamGP, nWorkers);

%% MC Plotting
% Convergence Plotting
plotStateMoments(const, times, xStats_MC, alphaStats_MC, 'MC');

% Plotting State Convergence At End time
tToPlotConvergenceIdx = find(times == const.tEnd);
plotConvergenceOfStateMoments(const, tToPlotConvergenceIdx, xStats_MC,...
                                alphaStats_MC, Nmc);

if fPlotSimData
    plotSimulationResults(rawData_MC.simData,[]);
end

%% Check Bi-Modality
plotMC_HistogramsInteractive(const, rawData_MC, nSamplesPerBin);

%% PCE Approximation
[Cmatrix, pceData] = FitPCE_ToMC(rawData_MC, maxTotalOrderPCE, times, fGP, dInputs, ...
                        stochasticParam, stochasticParamGP, const);

%% SOBOL Indices Generations
[tau] = CalcSobolIdx(Cmatrix, pceData, dInputs);

figure('Name','SOBOL Indices');
for iState = 1:const.nStates
    subplot(2,3,iState);
    hold on;
    for inputK = 1:dInputs
        plot(times, squeeze(tau(inputK,:,iState)), ...
            'DisplayName',pceData.variableNames{inputK});
    end
    grid on;
    xlabel('Time (s)');
    ylabel(const.stateNames{iState});
end
legend();
sgtitle('PCE SOBOL Indices vs Time For Each State');
