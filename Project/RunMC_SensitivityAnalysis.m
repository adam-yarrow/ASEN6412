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
Nmc = 50000; 
nSamplesPerBin = 100; % For plotting interactive histograms

% PCE Parameters
maxTotalOrderPCE = 7; % (p)
nWorkers = 6; % normally 6

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
                        stochasticParam, stochasticParamGP, const, nWorkers);

RMSE_Mean = mean(pceData.RMSE, 2);
RMSE_Sigma = var(pceData.RMSE, 0, 2).^0.5;

figure('Name',sprintf('Ci vs Index for End Time of Each State - p = %i', maxTotalOrderPCE));
for iState = 1:5
    subplot(2,3,iState)
    semilogy(abs(Cmatrix(:,end,iState)),'.');
    grid on;
    xlabel('Index')
    ylabel('|c_i|');
    title(const.stateNames{iState});
end
sgtitle(sprintf('Ci vs Index for End Time of Each State - p = %i', maxTotalOrderPCE))

% %% TODO - maybe compare PCE mean and variance estimates with MC response?
% %% TODO - try different levels of approximation and see what the response is like in terms of SOBOL and MSE
% %% TODO - generate larger MC samples proportional to higher order terms

%% SOBOL Indices Generations
[tau] = CalcSobolIdx(Cmatrix, pceData, dInputs);

figure('Name',sprintf('PCE SOBOL Indices vs Time For Each State - p = %i',maxTotalOrderPCE));
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
sgtitle(sprintf('PCE SOBOL Indices vs Time For Each State - p = %i',maxTotalOrderPCE));

