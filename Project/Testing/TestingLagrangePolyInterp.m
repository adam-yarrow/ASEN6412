% [y, w] = spquad(3,3);
% 
% nTimes = 10;
% nStates = 2;
% N = numel(w);
% 
% yInterp = y;
% fVals = randn(nTimes, nStates, N);
% 
% fInterp = LagrangeInterpND(yInterp, fVals, y);
% 
% fError = norm(fInterp(:) - fVals(:)) / norm(fVals(:));

%% Assuming ran RunForwardUQ.m for now (will move into their too)

% normalize input variables
[yNormMC, ~, ~] = GetNormalizedStdVariables(rawData_MC,...
                        stochasticParam);
[yNormGP_MC] = GetNormalizedGPVariables(rawData_MC, stochasticParamGP);
yNormMC = vertcat(yNormMC, yNormGP_MC); yNormMC = yNormMC';

% SC normalized
[yNormSC, ~, ~] = GetNormalizedStdVariables(rawData_SC,...
                        stochasticParam);
[yNormGP_SC] = GetNormalizedGPVariables(rawData_SC, stochasticParamGP);
yNormSC = vertcat(yNormSC, yNormGP_SC'); yNormSC = yNormSC';

% GOAL: If can recreate MC runs with lagrange interp accurately then dodgy
% hack is working
%% interp

xInterpMCfromSC = LagrangeInterpND(yNormMC, rawData_SC.x, yNormSC);

% errors
xError = sqrt(sum((xInterpMCfromSC - rawData_MC.x).^2, 3));

figure();
idxToPlot = 100;
for i = 1:5
    subplot(2,3,i)
    hold on;
    % plot(times, xError(i,:));
    plot(times, rawData_MC.x(i,:,idxToPlot),'DisplayName','Truth MC');
    plot(times, xInterpMCfromSC(i,:,idxToPlot),'DisplayName','Interped');
    grid on;
    xlabel('Time (s)');
    ylabel(const.stateNames{i});
end
legend();

figure();
for i = 1:5
    subplot(2,3,i)
    hold on;
    plot(times, xError(i,:));
    % plot(times, rawData_MC.x(i,:,idxToPlot),'DisplayName','Truth MC');
    % plot(times, xInterpMCfromSC(i,:,idxToPlot),'DisplayName','Interped');
    grid on;
    xlabel('Time (s)');
    ylabel(const.stateNames{i});
end
legend();