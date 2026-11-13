function [fTheta_Y, rawData] = DiscreteBayesFilter(theta, thetaPath, y,...
                                            fThetaPrior, ...
                                            controller, const, ...
                                            nWorkers, integrationWeights)
    
    %{
        N-dimensional discrete Bayes filter. Calculates a discrete
        probability density for the parameters of interest, given a set of
        measurements.

        Inputs:
            theta = N x d array of parameters to test. Assumes a tensorized
                grid has been made (or a sparse grid) and flattened. E.g. N =
                Nd^d, d = number of parameters
            thetaPaths = dx1 cell array defining the paths the theta
                parameters point to
            y = measurements nMeas x nTimes
            fThetaPrior = N x 1 array containing the prior on the theta
                parameters (again flattened)
            controller = controller function
            const = baseline constants
            integrationWeights = weights associated with each theta for
                integrating out theta to get fY(y) for the bayes rule
                denominator
        
        Outputs:
            fTheta_Y = posterior of size N x 1, corresponding to theta
            rawData = debugging structure
    %}
    fNormalizeAsProbability = false;
    if nargin < 8
        fNormalizeAsProbability = true;
    end
    if nargin < 7
        nWorkers = 0; % e.g. linear not parallel
    end
    
    %% Setup Storage
    [N, d] = size(theta);
    [nMeas, nTime] = size(y);

    fTheta_Y = zeros(N,1);    

    likelihoodVsTimeAndTheta = NaN(nTime, N);
    measInnov = NaN(nMeas, nTime, N);

    % tStart = tic;
    % wb = progressBar(0, N, tStart, [], 'Discrete Bayes Filter');
    WaitMessage = parfor_wait(N, 'Waitbar', true);
    parfor (iTheta = 1:N, nWorkers)
        currentTheta = theta(iTheta,:);
        [logLikelihood, ...
            likelihoodVsTimeAndTheta(:,iTheta),...
            measInnov(:,:,iTheta)] = logLikelihoodFcn(y, currentTheta, ...
                                            thetaPath, controller, const);

        % Update fTheta_Y vector
        fTheta_Y(iTheta) = logLikelihood + log(fThetaPrior(iTheta));

        % wb = progressBar(iTheta, N, tStart, wb, 'Discrete Bayes Filter');
        WaitMessage.Send;
    end

    WaitMessage.Destroy
    % Convert From Log to Normal Domain
    % fTot = LogSumExp(fTheta_Y);
    % fTheta_Y = exp(fTheta_Y - fTot);

    maxFTheta_Y = max(fTheta_Y);
    fTheta_Y = fTheta_Y - maxFTheta_Y; % rescale to avoid numerical issues
    fTheta_Y = exp(fTheta_Y);

    if fNormalizeAsProbability
        fTheta_Y = fTheta_Y ./ sum(fTheta_Y);
    else
        fY = sum(fTheta_Y .* integrationWeights); % Basic quadrature rule
        fTheta_Y = fTheta_Y ./ fY;
    end
    
    %% Save Raw Data
    rawData = struct();
    rawData.likelihoodVsTimeAndTheta = likelihoodVsTimeAndTheta;
    rawData.measInnov = measInnov;
end