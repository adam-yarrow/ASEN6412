function plotMCvsSC_StateMoments(const, times, stateMeanRelDiff, stateCovRelDiff, ...
                                    alphaMeanRelDiff, alphaVarRelDiff)
    
    figure('Name',sprintf('SC State Statistics Relative to MC Statistics Vs Time'));
    hStates = [];
    for iState = 1:const.nStates  
        hStates(end+1) = subplot(2,3,iState);
        hold on;
        stateSF = const.statePlottingSF(iState);
        plot(times(2:end), stateSF*stateMeanRelDiff(iState,2:end,end),'r',...
            'DisplayName','\mu');
        plot(times(2:end), stateSF*squeeze(stateCovRelDiff(iState,iState,2:end,end)),...
            'b','DisplayName','\sigma');
        
        legend();
        xlabel('Time (s)');
        ylabel(sprintf('%s Absolut Error (%s)', const.stateNames{iState}, const.statePlottingUnits{iState}));
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
    plot(times(2:end), stateSF*alphaMeanRelDiff(2:end),'r','DisplayName','\mu');
    plot(times(2:end), stateSF*alphaVarRelDiff(2:end),'b','DisplayName','\sigma');
                
    xlabel('Time (s)');
    ylabel(sprintf('%s Absolute Error (%s)', 'Alpha','deg'));
    grid on;   
    legend();
    ax = gca;
    ax.LineWidth = 2;  % Thicker axes
    ax.FontSize = 12;
    
    linkaxes(hStates,'x');
    sgtitle(sprintf('Difference Between SC and MC Statistics Vs Time'))

end
