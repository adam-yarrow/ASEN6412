function plotConvergenceOfStateMoments(const, tToPlotIdx, xStats, alphaStats, nMC)
  %{
    For plotting MC convergence of state mean and variance
  %}
    figure('Name',sprintf('State Moments Vs Nmc'));
    hStates = [];
    for iState = 1:const.nStates  
        hStates(end+1) = subplot(2,3,iState);
       
        dataExceptEnd = squeeze(xStats.cov(iState,iState,tToPlotIdx,1:end-1));
        dataEnd = xStats.cov(iState,iState,tToPlotIdx,end);
        dataRelDecrement = log10(abs((dataExceptEnd - dataEnd) ./ dataEnd));
        semilogx(1:nMC-1, dataRelDecrement,'DisplayName','Variance');

        hold on;
        dataExceptEnd = squeeze(xStats.mean(iState,tToPlotIdx,1:end-1));
        dataEnd = xStats.mean(iState,tToPlotIdx,end);
        dataRelDecrement = log10(abs((dataExceptEnd - dataEnd) ./ dataEnd));
        semilogx(1:nMC-1, dataRelDecrement,'DisplayName','Mean');
        
        xlabel('Number of Simulations (s)');
        ylabel(sprintf('log_{10}(|%s - %s_{ref}|/%s_{ref})', ...
            const.stateNames{iState},const.stateNames{iState},const.stateNames{iState}));
        grid on;  
        ax = gca;
        ax.LineWidth = 2;  % Thicker axes
        ax.FontSize = 12;
        legend();
    end
    
    % Alpha
    hStates(end+1) = subplot(2,3,6);
   
    dataExceptEnd = squeeze(alphaStats.var(tToPlotIdx,1:end-1));
    dataEnd = alphaStats.var(tToPlotIdx,end);
    dataRelDecrement = log10(abs((dataExceptEnd - dataEnd) ./ dataEnd));
    semilogx(1:nMC-1, dataRelDecrement,'DisplayName','Variance');
    
    hold on;
    dataExceptEnd = squeeze(alphaStats.mean(tToPlotIdx,1:end-1));
    dataEnd = alphaStats.mean(tToPlotIdx,end);     
    dataRelDecrement = log10(abs((dataExceptEnd - dataEnd) ./ dataEnd));
    semilogx(1:nMC-1, dataRelDecrement,'DisplayName','Mean');
        
    xlabel('Number of Simulations (s)');
    ylabel(sprintf('log_{10}(|%s - %s_{ref}|/%s_{ref})', ...
        '\alpha','\alpha','\alpha'));
    grid on;   
    ax = gca;
    ax.LineWidth = 2;  % Thicker axes
    ax.FontSize = 12;
    legend();
    
    linkaxes(hStates,'x');
    sgtitle(sprintf('Decrement in State Moments Relative to Final Solution vs Number of Samples'))



end
