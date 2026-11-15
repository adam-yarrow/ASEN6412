function [P] = LegendrePoly(y, order)
    %{
        Returns a normalized legendre polynomial evaluated at points y.
    %}

    % I have some conflicting paths to deal with... 
    w = warning('off', 'MATLAB:class:AtFolderPrecedence');     
    P = legendreP(order, y) / sqrt(2/(2*order + 1));
    warning(w);
    
end