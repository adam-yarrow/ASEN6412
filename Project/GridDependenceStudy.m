%{
    TODO - run different dT and show getting good accuracy 
     just running baseline case
%}

Startup();
const = ModelParams(); % Baseline Parameters
controller =  @(t,x,const) AlphaController(t, x, const, @(t) getReferenceAlpha(t));

dToptions = [0.1, 0.05, 0.01, 0.001, 0.0001];
storageTimes = 0:dToptions(1):const.tEnd;
nTimes = numel(storageTimes);

nDT = numel(dToptions);

xStorage = NaN(const.nStates, nTimes, nDT-1);
alphaStorage = NaN(nTimes, nDT-1);
simDataStorage = cell(nDT,1);

refX = [];
refAlpha = [];

dtNames = cellfun(@num2str, num2cell(dToptions), 'UniformOutput', false);

for iDT = 1:nDT
    const.dT = dToptions(iDT);
    simData = Simulation(controller, const);

    simDataStorage{iDT} = simData;

    [~, validIdx] = intersect(simData.times, storageTimes);

    if iDT ~= nDT
        xStorage(:,:,iDT) = simData.x(:,validIdx);
        alphaStorage(:,iDT) = calcAlpha( simData.x(:,validIdx));  
    else
        refX = simData.x(:,validIdx);
        refAlpha = calcAlpha( simData.x(:,validIdx));  
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
    for iDT = 1:nDT - 1
        plot(storageTimes, squeeze(xStorage(iState,:,iDT))*const.statePlottingSF(iState))
    
        xlabel('Time (s)');
        ylabel(sprintf('%s (%s)', const.stateName{iState}, const.statePlottingUnits{iState}));
        grid on;  
        legend();
        ax = gca;
        ax.LineWidth = 2;  % Thicker axes
        ax.FontSize = 12;
    end
end


% plotSimulationResults(simDataStorage, dtNames);

