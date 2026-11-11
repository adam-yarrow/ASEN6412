function [xStats, alphaStats, xStatsLogDiff, alphaStatsLogDiff] = ...
    CheckStochasticCollocConvergence(cvgLevels, ...
                Yref, refSimulationData, xStatsRef, alphaStatsRef,...
                const, controller, stochasticParam, stochasticParamGP,...
                comparisonTol)
    %{
        Checks the convergence of stochastic collocation as a function of
        level parameter. Assumes you have ran a reference simulation (Yref,
        refSimulationData) that is of a higher level parameter than the
        levels in cvgLevels.

        Inputs:
            cvgLevels = Nx1 array with the levels to check convergence on.
            Yref = Nxd vector of reference smolyak absisscas.
            refSimulationData = cell array containing the simData.
            xStatsRef = reference statistics
            alphaStatsRef = reference alphastats
            evaluations for the reference simulation.
            const = same constant parameters as reference.
            controller = same controller func as reference.
            stochasticParam = same stochastic param as reference.
            stochasticParamGP = for Gaussian Process of CMalpha

        Outputs:    
            xStats = structure containing mean and covariance stats for
            each case
            alphaStats = like xStats.
            
            
    %}
    % General Sizes
    [~, dInputs] = size(Yref);
    nTimes = numel(refSimulationData{1}.times);
    nLevels = numel(cvgLevels);

    % Storage
    xStats.mean = NaN(const.nStates, nTimes, nLevels);
    xStats.cov = NaN(const.nStates, const.nStates, nTimes, nLevels);
    xStats.corr = NaN(const.nStates, const.nStates, nTimes, nLevels);

    alphaStats.mean = NaN(nTimes, nLevels);
    alphaStats.var = NaN(nTimes, nLevels);

    xStatsLogDiff = xStats;
    alphaStatsLogDiff = alphaStats;

    % Run SC for each sub level
    for iLevel = 1:numel(cvgLevels)
        % Get current Y and W
        [Y, W] = spquad(dInputs, cvgLevels(iLevel));       

        % Find indices corresponding to the current level from the
        % reference
        [yInRefIdx, yInRefLoc] = ismembertol(Yref, Y, comparisonTol,'ByRows',true);
    
        % This removes all the zeros to get the ordering of how Yi appear 
        % in the refernce
        yInRefLoc(yInRefLoc == 0) = [];         
        [~, sortedIdx] = sort(yInRefLoc);


        % Pull out relevant sim data and sort according to yInRefLoc
        % (Probably dont need to sort cause normally in order anyway, but
        % safer this way)
        simData = refSimulationData(yInRefIdx);
        simData = simData(sortedIdx); % sorting

        % Run SC
        [xStats_i, alphaStats_i, ~] = StochasticColloc(const, Y, W, ...
                         controller, stochasticParam, stochasticParamGP, ...
                         simData);

        % Save Data
        xStats.mean(:,:,iLevel)  = xStats_i.mean;
        xStats.cov(:,:,:,iLevel)  = xStats_i.cov;
        xStats.corr(:,:,:,iLevel)  = xStats_i.corr;

        alphaStats.mean(:,iLevel)  = alphaStats_i.mean;
        alphaStats.var(:,iLevel)  = alphaStats_i.var;       

        % Generate Abs Log Differences Relative to Reference Traj
        xStatsLogDiff.mean = log10(abs(xStats.mean - xStatsRef.mean) ./ xStatsRef.mean);
        xStatsLogDiff.cov = log10(abs(xStats.cov - xStatsRef.cov) ./ xStatsRef.cov); % TODO - might want to do this on std dev?
        xStatsLogDiff.corr = log10(abs(xStats.corr - xStatsRef.corr) ./ xStatsRef.corr);

        alphaStatsLogDiff.mean = log10(abs(alphaStats.mean - alphaStatsRef.mean) ./ alphaStatsRef.mean);
        alphaStatsLogDiff.var = log10(abs(alphaStats.var - alphaStatsRef.var) ./ alphaStatsRef.var);

    end
end