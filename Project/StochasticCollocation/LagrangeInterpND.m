function fInterp = LagrangeInterpND(yInterp, fVals, y)
%{
    Interpolates the SC results at a specified y vector

    yInterp = M x d
    fVals = nStates x nTimes x N
    y = N x d
    
%}
    [M, d] = size(yInterp);
    [nStates, nTimes, N] = size(fVals);

    fInterp = zeros(nStates, nTimes, M);
    for iSample = 1:N
        lagrangeProduct = ones(M,1);
        for iDim = 1:d
            currentDimProd = LagrangePoly1D(y(:,iDim), ...
                                yInterp(:,iDim), iSample);
            lagrangeProduct = lagrangeProduct .* currentDimProd;
        end
               
        lagrangeProduct = reshape(lagrangeProduct, 1, 1, M);
        fInterp = fInterp + fVals(:,:,iSample) .* lagrangeProduct;
    end
end


