function [muU, varU] = CumulativeStats(u, Nmc)
    % Running Mean
    muU = cumsum(u,2) ./ (1:Nmc);

    % % Running Var
    % uSquared = cumsum(u.^2,2);
    % NmeanU = (1:Nmc) .* muU.^2;
    % varU = (uSquared - NmeanU) ./ ((1:Nmc) - 1);
    % varU(:,1) = NaN;
    % warning('cumulative variance might have issues with small variances relative to mean')


    % Welfords Update - should address my algorithm's issues
    c = u(:,1);
    d = u - c;
    S1 = cumsum(d,2);
    S2 = cumsum(d.^2,2);
    varU = (S2 - S1.^2 ./ Nmc) ./ (Nmc - 1);
    varU(:,1) = NaN;
end