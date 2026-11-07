function plotSimulationResults(simDataCell,simNames, refTrajFunc)
%{
    Plots the simulation response.
    simData can be a cell array
%}

const = ModelParams();

if nargin < 3
    refTrajFunc = [];
end

%% Ref Traj
if ~isempty(refTrajFunc)
    xr = refTrajFunc(simDataCell{1}.times);    
else
    xr = [];
end

nTraj = length(simDataCell);


%% Dynamic pressure
figure('Name',sprintf('Dynamic Pressure'));
hold on;
for iTraj = 1:nTraj
    simData = simDataCell{iTraj};
    if ~isempty(simNames)
        plot(simData.times, simData.qBar/1000,'DisplayName',simNames{iTraj},'LineWidth',1.5);
    else
        plot(simData.times, simData.qBar/1000,'LineWidth',1.5);
    end
end
if ~isempty(simNames)
    legend();
end
xlabel('Time (s)');
ylabel('Dynamic Pressure (kPa)');
grid on;
title(sprintf('Qbar vs Time'));
ax = gca;
ax.LineWidth = 2;  % Thicker axes
ax.FontSize = 12;

%% Control Inputs
figure('Name',sprintf('Control inputs'));
hold on;
iInput = 1;
plotState(simDataCell, 'u', iInput, simNames, const.inputNames{iInput},...
            const.inputPlottingUnits{iInput}, const.inputPlottingSF(iInput),[]);
grid on;   
ax = gca;
ax.LineWidth = 2;  % Thicker axes
ax.FontSize = 12;
title(sprintf('Inputs vs Time'));

%% States
figure('Name',sprintf('States'));
hStates = [];
for iState = 1:const.nStates
    if ~isempty(xr)
        xrCurrent = xr(iState,:);
    else
        xrCurrent = [];
    end

    hStates(end+1) = subplot(2,3,iState);
    hold on;
    plotState(simDataCell, 'x', iState, simNames, ...
        const.stateNames{iState}, const.statePlottingUnits{iState},...
        const.statePlottingSF(iState), xrCurrent);
    % legend();
end

% Alpha
hStates(end+1) = subplot(2,3,6);
hold on;
for iTraj = 1:numel(simDataCell)
    simData = simDataCell{iTraj};
    data = calcAlpha(simData.x);
    t = simData.times;
    plot(t,data*rad2deg(1),'LineWidth',1.5);%,'DisplayName',simNames{iTraj});
    hold on;
end
xlabel('Time (s)');
ylabel(sprintf('%s (%s)', 'Alpha', 'deg'));
grid on;   
if ~isempty(simNames)
    legend(simNames{:});
end
ax = gca;
ax.LineWidth = 2;  % Thicker axes
ax.FontSize = 12;

linkaxes(hStates,'x');
sgtitle(sprintf('States vs Time'))


end


function plotState(simDataCell, dataName, stateIdx, simNames, stateName, stateUnits,SF, xr)
    for iTraj = 1:numel(simDataCell)
        simData = simDataCell{iTraj};
        data = simData.(dataName)(stateIdx,:);
        t = simData.times;
        if ~isempty(simNames)
            plot(t,data*SF,'DisplayName',simNames{iTraj},'LineWidth',1.5);
        else
            plot(t,data*SF,'LineWidth',1.5);
        end
    end

    if ~isempty(xr)
        plot(t, xr*SF,'r--','DisplayName','Reference Trajectory','LineWidth',1.5);
    end

    xlabel('Time (s)');
    ylabel(sprintf('%s (%s)', stateName, stateUnits));
    grid on;  
    if ~isempty(simNames)
        legend();
    end
    ax = gca;
    ax.LineWidth = 2;  % Thicker axes
    ax.FontSize = 12;
end