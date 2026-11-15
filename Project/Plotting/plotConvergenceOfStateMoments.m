function plotConvergenceOfStateMoments(const, tToPlotIdx, xStats, alphaStats, nMC)
  %{
    For plotting MC convergence of state mean and variance
  %}
%% Mean
plotShit(const, tToPlotIdx, xStats, alphaStats, nMC, 'Mean');
%% Variance
plotShit(const, tToPlotIdx, xStats, alphaStats, nMC, 'Variance');

end

%% Helper
function plotShit(const, tToPlotIdx, xStats, alphaStats, nMC, overallName)
    figure('Name',sprintf('State %s Vs Nmc',overallName));
    hStates = [];
    for iState = 1:const.nStates  
        hStates(end+1) = subplot(2,3,iState);
        % hold on;
        if strcmpi(overallName, 'Variance')
            dataExceptEnd = squeeze(xStats.cov(iState,iState,tToPlotIdx,1:end-1));
            dataEnd = xStats.cov(iState,iState,tToPlotIdx,end);
        else
            dataExceptEnd = squeeze(xStats.mean(iState,tToPlotIdx,1:end-1));
            dataEnd = xStats.mean(iState,tToPlotIdx,end);
        end        
        dataRelDecrement = log10(abs((dataExceptEnd - dataEnd) ./ dataEnd));
        semilogx(1:nMC-1, dataRelDecrement);
        
        xlabel('Number of Simulations (s)');
        ylabel(sprintf('log_{10}(|%s - %s_{ref}|/%s_{ref})', ...
            const.stateNames{iState},const.stateNames{iState},const.stateNames{iState}));
        grid on;  
        ax = gca;
        ax.LineWidth = 2;  % Thicker axes
        ax.FontSize = 12;
    end
    
    % Alpha
    hStates(end+1) = subplot(2,3,6);
    if strcmpi(overallName, 'Variance')
        dataExceptEnd = squeeze(alphaStats.var(tToPlotIdx,1:end-1));
    dataEnd = alphaStats.var(tToPlotIdx,end);
    else
        dataExceptEnd = squeeze(alphaStats.mean(tToPlotIdx,1:end-1));
        dataEnd = alphaStats.mean(tToPlotIdx,end);
    end   
    dataRelDecrement = log10(abs((dataExceptEnd - dataEnd) ./ dataEnd));
    semilogx(1:nMC-1, dataRelDecrement);
        
    xlabel('Number of Simulations (s)');
    ylabel(sprintf('log_{10}(|%s - %s_{ref}|/%s_{ref})', ...
        '\alpha','\alpha','\alpha'));
    grid on;   
    ax = gca;
    ax.LineWidth = 2;  % Thicker axes
    ax.FontSize = 12;
    
    linkaxes(hStates,'x');
    sgtitle(sprintf('Decrement in State %s Relative to Final Solution vs Number of Samples',overallName))


end