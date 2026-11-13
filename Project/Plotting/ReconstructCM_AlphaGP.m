function CM_alpha = ReconstructCM_AlphaGP(Y, times, dataKLE, idxToPlot)
    %{
        Plots the CM_alpha values as a function of time for different MC
        instances to see how they behave.
    %}

    if nargin < 4
        fPlotAll = true; 
    else 
        fPlotAll = false; 
    end
    
    phi = MakePhiMatrix(times, dataKLE);    
    
    CM_alpha = getRV_KLE_Fast(times, dataKLE, phi, Y);

    figure('Name','CM_Alpha Gaussian Process');
    if fPlotAll
        hold on;
        for i = 1:size(CM_alpha,2)    
            plot(times, CM_alpha(:,i));
        end
    else
        plot(times, CM_alpha(:,idxToPlot));
    end
    grid on;
    xlabel('Time (s)');
    ylabel('C_M_{\alpha}');
    title('C_M_{\alpha} vs Time')
    

    
end
