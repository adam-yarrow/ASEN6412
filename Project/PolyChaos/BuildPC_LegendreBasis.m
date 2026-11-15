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

    % Pre-Compute univariate basis functions Once
    univariateBasisEvals = NaN(p+1, d);
    for k = 0:p
        %% TODO - might make this universal and provide option to specify poly type here in future
        %% Would need to fuck with pre-normalization though
        univariateBasisEvals(k+1,:) = LegendrePoly(y,k);
    end
    
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