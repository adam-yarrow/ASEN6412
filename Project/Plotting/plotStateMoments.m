function plotStateMoments(const, times, xStats, alphaStats, dataSrc)
    
    figure('Name',sprintf('State Statistics Vs Time - %s', dataSrc));
    hStates = [];
    for iState = 1:const.nStates  
        hStates(end+1) = subplot(2,3,iState);
        hold on;
        stateSF = const.statePlottingSF(iState);
        sigma = squeeze(xStats.cov(iState,iState,2:end,end)).^0.5;
        plot(times(2:end), stateSF*xStats.mean(iState,2:end,end),'k','DisplayName','\mu');
        plot(times(2:end), stateSF*(squeeze(xStats.mean(iState,2:end,end)) - sigma'),'r--','DisplayName','\mu - 1\sigma');
        plot(times(2:end), stateSF*(squeeze(xStats.mean(iState,2:end,end)) + sigma'),'b--','DisplayName','\mu + 1\sigma');
        
        legend();
        xlabel('Time (s)');
        ylabel(sprintf('%s (%s)', const.stateNames{iState}, const.statePlottingUnits{iState}));
        grid on;  
        legend();
        ax = gca;
        ax.LineWidth = 2;  % Thicker axes
        ax.FontSize = 12;
    end
    
    % Alpha
    hStates(end+1) = subplot(2,3,6);
    hold on;
    stateSF = rad2deg(1);
    sigma = alphaStats.var(2:end,end).^0.5;
    plot(times(2:end), stateSF*alphaStats.mean(2:end,end),'k','DisplayName','\mu');
    plot(times(2:end), stateSF*(alphaStats.mean(2:end,end) - sigma),'r--','DisplayName','\mu - 1\sigma');
    plot(times(2:end), stateSF*(alphaStats.mean(2:end,end) + sigma),'b--','DisplayName','\mu + 1\sigma');
        
    xlabel('Time (s)');
    ylabel(sprintf('%s (%s)', 'Alpha', 'deg'));
    grid on;   
    legend();
    ax = gca;
    ax.LineWidth = 2;  % Thicker axes
    ax.FontSize = 12;
    
    linkaxes(hStates,'x');
    sgtitle(sprintf('State Statistics vs Time - %s', dataSrc))

end
