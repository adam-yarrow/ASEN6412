function [phiVsP] = BuildPC_LegendreBasis(y, validBasisCombos)
    %{
        Builds and evaluates the polynomial basis functions for a given
        vector of random values (y). To be used for finding PCE
        coefficients for SOBOL indice generation.

        y = d x 1
        validBasisCombos = (P+1) x d array where each row is a
        different basis combination (e.g. [j1, j2, ... jd])

        phiVsP = Row vector of size 1xP+1, P+1 = (p+d)!/(p!d!) containing
        the basis function (up to total order p) evaluated at the y values.
    %}
    [nBasisFunc, d] = size(validBasisCombos);
    p = max(validBasisCombos,[],'all');

    % Pre-Compute univariate basis functions once using 3 term recurrence
    univariateBasisEvals = NaN(p+1, d);
    normConstant = NaN(p+1,1);
    y = y(:)'; % forcing row vector

    for k = 0:p
        normConstant(k+1) = sqrt((2*k + 1));  % OR is it sqrt((2*k + 1)/2) - this is for no prob density i think?
        if k == 0
            univariateBasisEvals(k+1,:) = ones(1,d);
        elseif k == 1
            univariateBasisEvals(k+1,:) = y;    
        else
            i = k - 1;
            % phi_i = univariateBasisEvals(i+1,:);
            % phi_it1 = univariateBasisEvals(i-1 + 1,:);
            % solving for phi_ip1
            univariateBasisEvals(k+1,:) = 1/(i+1) * ...
                    ((2*i + 1)*y .* univariateBasisEvals(i+1,:) ...
                        - i*univariateBasisEvals(i-1 + 1,:));
        end
    end
    % Normalize
    univariateBasisEvals = univariateBasisEvals .* normConstant;

    % Evaluate Basis Functions
    phiVsP = ones(1, nBasisFunc);
    % TODO - determine if i can re-use anything here?? OR how to avoid
    % loops
    for iBasis = 1:nBasisFunc
        for iDim = 1:d
            jk = validBasisCombos(iBasis,iDim);
            phiVsP(1, iBasis) = phiVsP(1, iBasis) * univariateBasisEvals(jk+1, iDim);
        end
    end    
end