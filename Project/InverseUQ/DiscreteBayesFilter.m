function [fTheta_Y, rawData] = DiscreteBayesFilter(theta, thetaPath, y,...
                                            fThetaPrior, ...
                                            controller, const)
    
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
        
        Outputs:
            fTheta_Y = posterior of size N x 1, corresponding to theta
            rawData = debugging structure
    %}
%% NOTE: Could maybe try a sparse grid on this and then fit a response surface to it?

    %% Setup Storage
    [N, d] = size(theta);
    [nMeas, nTime] = size(y);

    fTheta_Y = zeros(N,1);    

    likelihoodVsTimeAndTheta = NaN(nTime, N);
    measInnov = NaN(nMeas, nTime, N);

    tStart = tic;
    wb = progressBar(0, N, tStart, [], 'Discrete Bayes Filter');
    for iTheta = 1:N
        currentTheta = {theta(iTheta,:)};
        [logLikelihood, ...
            likelihoodVsTimeAndTheta(:,iTheta),...
            measInnov(:,:,iTheta)] = logLikelihoodFcn(y, currentTheta, ...
                                            thetaPath, controller, const);

        % Update fTheta_Y vector
        fTheta_Y(iTheta) = logLikelihood + log(fThetaPrior(iTheta));

        wb = progressBar(iTheta, N, tStart, wb, 'Discrete Bayes Filter');
    end

    % Convert From Log to Normal Domain
    % fTot = LogSumExp(fTheta_Y);
    % fTheta_Y = exp(fTheta_Y - fTot);

    maxFTheta_Y = max(fTheta_Y);
    fTheta_Y = fTheta_Y - maxFTheta_Y; % rescale to avoid numerical issues
    fTheta_Y = exp(fTheta_Y);

    %% FIXME - temporarily normalizing to a probability (not a density)
    fTheta_Y = fTheta_Y ./ sum(fTheta_Y);

    % Normalize into a Proper Density (e.g. sums to 1.0) - TODO could use
    % LogSumExp here?
    % fTheta_Y = fTheta_Y ./ trapz(theta, fTheta_Y);

    warning('need to work out weighting here, if use sparse grids could use smolyak weights for integrating?')


    %% Save Raw Data
    rawData = struct();
    rawData.likelihoodVsTimeAndTheta = likelihoodVsTimeAndTheta;
    rawData.measInnov = measInnov;
end