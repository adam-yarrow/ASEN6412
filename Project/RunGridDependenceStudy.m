%{
    This script runs a mesh dependence study on my ODE solver to determine
    what relative error tolerance is required for an adequatly converged
    solution.

    Kind of a farce as by specifying a relative error tolerance that
    implicity fixes the precision, but this just demonstrates the impact
    the tolerance has on all of the states.
%}

Startup();
const = ModelParams(); % Baseline Parameters
controller =  @(t,x,const) AlphaController(t, x, const, @(t) getReferenceAlpha(t));

relErrorTol = [1E-1, 1E-2, 1E-3, 1E-4, 1E-5, 1E-6]; %[0.1, 0.05, 0.01, 0.001, 0.0001];
times = 0:const.dT:const.tEnd;
nTimes = numel(times);

nRelTol = numel(relErrorTol);

xStorage = NaN(const.nStates, nTimes, nRelTol-1);
alphaStorage = NaN(nTimes, nRelTol-1);
simDataStorage = cell(nRelTol,1);


refX = [];
refAlpha = [];

relTolNames = cell(nRelTol,1); % cellfun(@num2str, num2cell(relErrorTol), 'UniformOutput', false);

for i = 1:nRelTol

    const.odeOpts.relTol = relErrorTol(i);
    simData = Simulation(controller, const);
    simDataStorage{i} = simData;

    xInterp = interp1(simData.times, simData.x', times)';
    alphaInterp = interp1(simData.times, calcAlpha(simData.x)', times)';

    relTolNames{i} = sprintf('RelTol = 1E%i', log10(relErrorTol(i)));
    if i ~= nRelTol
        xStorage(:,:,i) = xInterp;
        alphaStorage(:,i) = alphaInterp;
    else
        refX = xInterp;
        refAlpha = alphaInterp;
    end
end

%% Plotting Convergence
% Calculate Log Decrement
logRelDiffX = log10(abs((xStorage - refX)./ refX));
logRelDiffAlpha = log10(abs((alphaStorage - refAlpha)./ refAlpha));

figure('Name','Mesh Dependence Study');
hStates = [];
for iState = 1:const.nStates
    hStates(end+1) = subplot(2,3,iState);
    hold on;
    for i = 1:nRelTol - 1
        plot(times, squeeze(logRelDiffX(iState,:,i)),'DisplayName',relTolNames{i});
    
        xlabel('Time (s)');
        ylabel(sprintf('%s', const.stateNames{iState}));
        grid on;  
        legend();
        ax = gca;
        ax.LineWidth = 2;  % Thicker axes
        ax.FontSize = 12;
    end
end

hStates(end+1) = subplot(2,3,6);
hold on;
for i = 1:nRelTol - 1
    plot(times, squeeze(logRelDiffAlpha(:,i)),'DisplayName',relTolNames{i});

    xlabel('Time (s)');
    ylabel(sprintf('%s', '\alpha'));
    grid on;  
    legend();
    ax = gca;
    ax.LineWidth = 2;  % Thicker axes
    ax.FontSize = 12;
end

titleText = sprintf(['Logarithm of the Relative Error Between Different ODE45', ...
    'Integration Error Tolerances with a Reference of 1E%i'], log10(relErrorTol(end)));
sgtitle(titleText);

plotSimulationResults(simDataStorage, relTolNames);

