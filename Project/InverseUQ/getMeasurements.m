function y = getMeasurements(const, x)
    %{
        Generates additive white gaussian noise for measurement model.
    %}
    C = const.measModel.C;
    Sxx = chol(const.measModel.Cxx,'lower');

    y = C*x + Sxx*randn(size(x));    
end