function T = BuildCmAlphaSensitivityMatrix(t, nParams)
    %{
        Polynomial sensitivity matrix
    %}
    T = zeros(numel(t), nParams);
    for i = 1:nParams
        T(:,i) = t.^(i-1);
    end
end