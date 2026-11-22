function plotCvgSC_Trajectories(const, times, xStats, alphaStats, fLogDiff, ...
            cvglevelNames, dataNameState, dataNameAlpha)

figure('Name',sprintf('Convergence in Trajectories - State %s', dataNameState));
hStates = [];
for iState = 1:const.nStates
    hStates(end+1) = subplot(2,3,iState);
    hold on;
    plotState(xStats, times, dataNameState, iState, cvglevelNames, ...
        const.stateNames{iState}, const.statePlottingUnits{iState},...
        const.statePlottingSF(iState), fLogDiff);
    % legend();
end

% Alpha
hStates(end+1) = subplot(2,3,6);
hold on;
if strcmpi(dataNameAlpha,'mean')
    SF = rad2deg(1);
else
    SF = rad2deg(1)^2;
end

if fLogDiff
    SF = 1;        
end
for iTraj = 1:numel(cvglevelNames)
    data = alphaStats.(dataNameAlpha)(:,iTraj);
    
    plot(times,data*SF,'DisplayName', ...
        sprintf('q = %s',cvglevelNames{iTraj}),'LineWidth',1.5);

end

xlabel('Time (s)');
if fLogDiff
    ylabel(sprintf('log_{10}(|(%s - %s_{ref})/%s_{ref}|)', '\alpha','\alpha', '\alpha'));
else
    if strcmpi(dataNameAlpha,'mean')
        ylabel(sprintf('%s (%s)', '\alpha', 'deg'));
    else
        ylabel(sprintf('%s (%s)', '\alpha', 'deg^2'));
    end
    
end

grid on;   
legend();
ax = gca;
ax.LineWidth = 2;  % Thicker axes
ax.FontSize = 12;

linkaxes(hStates,'x');
sgtitle(sprintf('Convergence in State %s vs Time and Level Parameter', dataNameState))

end


function plotState(statsArray, times, dataName, stateIdx, simNames, ...
    stateName, stateUnits, SF, fLogDiff)
        
    if strcmpi(dataName,'cov')
        SF = SF^2;
    end

    if fLogDiff
        SF = 1;
    end

    for iTraj = 1:numel(simNames)
        if strcmpi(dataName, 'mean')
            data = statsArray.(dataName)(stateIdx,:,iTraj);
        else
            data = squeeze(statsArray.(dataName)(stateIdx, stateIdx,:,iTraj));
        end
                       
        plot(times,data*SF,'DisplayName', ...
            sprintf('q = %s',simNames{iTraj}),'LineWidth',1.5);
 
    end

    xlabel('Time (s)');
    if fLogDiff
        ylabel(sprintf('log_{10}((|%s - %s_{ref})/%s_{ref}|)', stateName, stateName, stateName));
    else
        if strcmpi(dataName,'mean')
            ylabel(sprintf('%s - %s (%s)', stateName, dataName, stateUnits));
        else
            ylabel(sprintf('%s - %s (%s)^2', stateName, dataName, stateUnits));
        end
    end
    grid on;  
    % legend();
    ax = gca;
    ax.LineWidth = 2;  % Thicker axes
    ax.FontSize = 12;
end