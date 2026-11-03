function plotUxVsNpoints(titleName, x, muU, varU, Npoints, fUseLogPoints)
    %{
        This plots the mean and variance response of u(x) for different
        number of integration or MC points.
    %}
    figure('Name',sprintf('Mean and Var Response - %s', titleName));
    h1 = subplot(2,1,1);
    grid on;
    hold on;
    xlabel('x');
    ylabel('E[u(x)]');
    
    h2 = subplot(2,1,2);
    grid on;
    hold on;
    xlabel('x');
    ylabel('Var(u(x))');
    
    for i = 1:numel(Npoints)

        if fUseLogPoints
            currentName = sprintf('N = 1E%i',log10(Npoints(i)));
        else
            currentName = sprintf('N = %i',Npoints(i));
        end

        plot(h1, x, muU(:,i),'DisplayName',...
            currentName);
        plot(h2, x, varU(:,i),'DisplayName',...
            currentName);
    end
    
    sgtitle(sprintf('Mean and Variance of u(x) vs increasing %s Points', titleName))
    legend(h1);
    legend(h2);
end