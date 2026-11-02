function plotUxVsNpointsDifference(titleName, x, muU, varU, Npoints, fUseLogPoints)
    %{ 
        NOTE: treats the last entry (largest Npoints) as the reference
    %}

    figure('Name',sprintf('Difference - Mean and Var Response - %s', titleName));
    h1 = subplot(2,1,1);
    grid on;
    hold on;
    xlabel('x');
    ylabel('log_{10}(E[u(x)] - E[u_{ref}(x)])');
    
    h2 = subplot(2,1,2);
    grid on;
    hold on;
    xlabel('x');
    ylabel('log_{10}(Var(u(x)) - Var(u_{ref}(x)))');
    
    for i = 1:numel(Npoints) - 1

        if fUseLogPoints
            currentName = sprintf('N = 1E%i',log10(Npoints(i)));
        else
            currentName = sprintf('N = %i',Npoints(i));
        end

        plot(h1, x, log10(abs(muU(:,i) - muU(:,end))) ,'DisplayName',...
            currentName);
        plot(h2, x, log10(abs(varU(:,i) - varU(:,end))),'DisplayName',...
            currentName);
    end
    
    sgtitle(sprintf('Difference in Mean and Variance of u(x) vs increasing %s Points', titleName))
    legend(h1);
    legend(h2);
end