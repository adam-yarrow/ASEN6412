function [yNorm, dLocal, varNames] = GetNormalizedStdVariables(rawData_MC, ...
                                    stochasticParam)
    y = [rawData_MC.rawStochasticValues{:}]'; % d x N (N simulations)
    dLocal = size(y,1); % NOTE: only local cause GP handled separately
    varNames = cell(dLocal, 1);

    % Map y from [a,b]-> [-1,1] for PCE basis
    yNorm = NaN(size(y));
    for iDim = 1:dLocal
        b = stochasticParam{iDim}.max;
        a = stochasticParam{iDim}.min;
    
        yNorm(iDim,:) =  2*(y(iDim,:) - (b+a)/2)/(b-a);
        
        varNames{iDim} = stochasticParam{iDim}.path;
    end


end