function [xStats, alphaStats] = CalculateMC_StatisticsFast(xStorage, alphaStorage)
    %{
        Calculates the cumulative statistics for a MC simulation.
        mean = nStates x nTimes x nMC
        cov and corr = nStates x nStates x nTimes x nMC

        Partially AI Generated: Took original code and asked AI to improve
        efficiency.
    %}
    [nStates, nTimes, nMC] = size(xStorage);

    xStats.mean = NaN(nStates, nTimes, nMC);
    xStats.cov  = NaN(nStates, nStates, nTimes, nMC);
    xStats.corr = NaN(nStates, nStates, nTimes, nMC);
    alphaStats.mean = NaN(nTimes, nMC);
    alphaStats.var  = NaN(nTimes, nMC);

    % Running state (Welford)
    mu  = zeros(nStates, nTimes);
    M2  = zeros(nStates, nStates, nTimes);
    aMu = zeros(nTimes, 1);
    aM2 = zeros(nTimes, 1);

    diagIdx = 1:(nStates + 1):nStates^2;   % linear indices of the diagonal

    for i = 1:nMC
        % ---- states: update running mean and scatter matrix ----
        x     = xStorage(:,:,i);
        d     = x - mu;                    % deviation from old mean
        mu    = mu + d / i;
        d2    = x - mu;                    % deviation from new mean
        M2    = M2 + pagemtimes(permute(d,  [1 3 2]), 'none', ...
                                permute(d2, [1 3 2]), 'transpose');
        xStats.mean(:,:,i) = mu;

        % ---- alpha ----
        a    = alphaStorage(:,i);
        da   = a - aMu;
        aMu  = aMu + da / i;
        aM2  = aM2 + da .* (a - aMu);
        alphaStats.mean(:,i) = aMu;
        alphaStats.var(:,i)  = aM2 / max(i - 1, 1);   % var of 1 sample = 0, as in MATLAB

        % ---- cov / corr (undefined for a single sample) ----
        if i > 1
            xStats.cov(:,:,:,i) = M2 / (i - 1);

            v = reshape(M2, nStates^2, nTimes);
            v = reshape(v(diagIdx,:), nStates, 1, nTimes);   % variances (unscaled)
            xStats.corr(:,:,:,i) = M2 ./ sqrt(v .* permute(v, [2 1 3]));
        end
    end
end