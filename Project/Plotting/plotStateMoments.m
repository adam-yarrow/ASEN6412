function plotStateMoments(const, times, xStats, alphaStats)
    
    figure('Name',sprintf('State Statistics Vs Time'));
    hStates = [];
    for iState = 1:const.nStates  
        hStates(end+1) = subplot(2,3,iState);
        hold on;
        stateSF = const.statePlottingSF(iState);
        sigma = squeeze(xStats.cov(iState,iState,:,end)).^0.5;
        plot(times, stateSF*xStats.mean(iState,:,end),'k','DisplayName','\mu');
        plot(times, stateSF*(squeeze(xStats.mean(iState,:,end)) - sigma'),'r--','DisplayName','\mu - 1\sigma');
        plot(times, stateSF*(squeeze(xStats.mean(iState,:,end)) + sigma'),'b--','DisplayName','\mu + 1\sigma');
        
        legend();
        xlabel('Time (s)');
        ylabel(sprintf('%s (%s)', const.stateNames{iState}, const.stateUnits{iState}));
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
    sigma = alphaStats.var.^0.5;
    plot(times, stateSF*alphaStats.mean(:,end),'k','DisplayName','\mu');
    plot(times, stateSF*(alphaStats.mean(:,end) - sigma),'r--','DisplayName','\mu - 1\sigma');
    plot(times, stateSF*(alphaStats.mean(:,end) + sigma),'b--','DisplayName','\mu + 1\sigma');
        
    xlabel('Time (s)');
    ylabel(sprintf('%s (%s)', 'Alpha', 'deg'));
    grid on;   
    legend();
    ax = gca;
    ax.LineWidth = 2;  % Thicker axes
    ax.FontSize = 12;
    
    linkaxes(hStates,'x');
    sgtitle(sprintf('State Statistics vs Time'))

end
