function [tau] = CalcSobolIdx(Cmatrix, pceData, d)
    %{
        Calculates the SOBOL indices associated with the PCE approximation
        of the MC simulation.

        Only looks at the total sensitivity index of the solution with
        respect to a given input yk.

        Cmatrix = (P+1) x Ntimes x nStates
        pceData = struct from fitting PCE
        d = Number of inputs        
    %}
    [~, nTimes, nStates] = size(Cmatrix);
    jk = pceData.basisOrders_jk;

    % Tau is a function of time, input variable and state
    tau = NaN(d,nTimes, nStates);
   
    %% Get total variance
    % Total Var = nTimes x nStates
    totalVar = sum(Cmatrix(2:end,:,:).^2, 1);
    
    %% Calculate Total Sensitivity Index
    for k = 1:d
        % Sum over Ci^2 for any cases where i_k > 0
        basisIdxForSum = jk(:,k) > 0;
        varContribYk = sum(Cmatrix(basisIdxForSum).^2,1);
        tau(k,:,:) = varContribYk ./ totalVar;        
    end   
    
end