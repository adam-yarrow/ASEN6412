%{
    TODO - run different dT and show getting good accuracy 
     just running baseline case
%}

Startup();
const = ModelParams(); % Baseline Parameters
controller =  @(t,x,const) AlphaController(t, x, const, @(t) getReferenceAlpha(t));

dToptions = [0.1, 0.05, 0.01, 0.001];

nDT = numel(dToptions);


% xStorage = NaN(const.nStates, nTimes, nDT);
% alphaStorage = NaN(nTimes, nDT);
simDataStorage = cell(nDT,1);

for iDT = 1:numel(dToptions)
    const.dT = dToptions(iDT);
    simData = Simulation(controller, const);

    simDataStorage{iDT} = simData;

    % xStorage(:,:,iDT) = simData.x;
    % alphaStorage(:,iDT) = calcAlpha(simData.x);    
    % WOULD NEED TO INTERPOLATE TO KEEP ON CONSISTENT GRID
end

cvgLevelsSCnames = cellfun(@num2str, num2cell(dToptions), 'UniformOutput', false);

plotSimulationResults(simDataStorage, cvgLevelsSCnames);


% Plot Errors - TODO generate log decrease errors
