function plotMC_Histograms(const, rawData, timeToPlotAt)
    %{
        Plots histogram of state response at a given time
    %}

    timeIdxToPlotAt = find(rawData.simData{1}.times == timeToPlotAt);

    figure('Name',sprintf('MC State Variable Histograms at T = %2.2f',timeToPlotAt));

    % normal states
    for iState = 1:const.nStates
        subplot(2,3,iState);
        data = const.statePlottingSF(iState) * ...
            squeeze(rawData.x(iState,timeIdxToPlotAt,:));
        histogram(data,'Normalization','probability','NumBins',numel(data)/50); % scaling number of bins by number of samples
        grid on;
        xlabel(sprintf('%s (%s)', const.stateNames{iState}, const.statePlottingUnits{iState}));
        ylabel('Probability');
    end

    % Alpha
    subplot(2,3,6);
    data = rad2deg(1) * ...
        squeeze(rawData.alpha(timeIdxToPlotAt,:));
    histogram(data,'Normalization','probability','NumBins',numel(data)/50); % scaling number of bins by number of samples
    grid on;
    xlabel(sprintf('%s (%s)', 'Alpha', 'deg'));
    ylabel('Probability');

    sgtitle(sprintf('MC State Variable Histograms at T = %2.2f', timeToPlotAt));
end



%{
Should try to make this work but very finiky

iState = 1; % state 2 - 0.0005

%% TODO - find a way to make this into a nicer surface?

stateMin = min(rawData_MC.x(iState,:,:),[],'all');
stateMax = max(rawData_MC.x(iState,:,:),[],'all');
nTimes = numel(times);

nXpoints = 1000;
xGrid = linspace(stateMin, stateMax, nXpoints);
fXest = NaN(nTimes,nXpoints);

for iTime = 301:nTimes
    [fXest(iTime, :)] = ksdensity(squeeze(rawData_MC.x(iState,iTime,:)), ...
        xGrid,'Bandwidth',1);
end

[X, T] = meshgrid(xGrid, times);

figure();
surf(X, T, fXest,'EdgeColor','none');
grid on;
xlabel('x')
ylabel('Time (s)');
zlabel('f_X(x)');

%%
%}