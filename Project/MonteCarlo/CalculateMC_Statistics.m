function [xStats, alphaStats] = CalculateMC_Statistics(xStorage, alphaStorage)
    %{
        Calculates the cumulative statistics for a MC simulation.
        
        mean = nStates x nTimes x nMC
        cov and corr = nStates x nStates x nTimes x nMC
    %}
    [nStates,nTimes,nMC] = size(xStorage);

    xStats.mean = NaN(size(xStorage));
    xStats.cov = NaN(nStates, nStates, nTimes, nMC);
    xStats.corr = NaN(nStates, nStates, nTimes, nMC);
            
    alphaStats.mean = NaN(nTimes, nMC);
    alphaStats.var = NaN(nTimes, nMC);

    % TODO - this isnt efficient, consider changing to a cumulative
    % approach (which will be a pain for covariance)
    for i = 1:nMC
        localX = xStorage(:,:,1:i);
        xStats.mean(:,:,i) = mean(localX, 3);
        
        dx = localX - xStats.mean(:,:,i);
        dxp = permute(dx, [1 3 2]); % nStates x nMC x nTimes
        xStats.cov(:,:,:,i) = pagemtimes(dxp, 'none', dxp, 'transpose') / (i - 1);

        for k = 1:nTimes
            xStats.corr(:,:,k,i) = corrcoef(reshape(localX(:,k,:), nStates, i).'); 
        end

        localAlpha = alphaStorage(:,1:i);
        alphaStats.mean = mean(localAlpha,2);
        alphaStats.var = var(localAlpha,0,2);
        
    end

end