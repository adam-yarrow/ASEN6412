function [validBasisCombos] = GetValidPC_BasisOptions(d,p)
    %{
        Generates the list of valid PCE basis options (k values) for a given
        dimensionality and maximum total order.
    %}

     % first build all possible combinations of y with individual order p
    Yoptions = cell(d,1);
    for iBasis = 1:d
        Yoptions{iBasis} = 0:p;
    end

    allBasisCombos = table2array(combinations(Yoptions{:}));

    % Now screen only options that satisfy our total order requirement
    validBasisCombosIdx = sum(allBasisCombos,2) <= p; % 1-norm along rows
    validBasisCombos = allBasisCombos(validBasisCombosIdx,:);
end