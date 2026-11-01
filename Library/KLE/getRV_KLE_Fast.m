function [x] = getRV_KLE_Fast(t, dataKLE, phiMatrix, Y)
 %{
    Returns the random value of the stochastic process at
    location/time t.

    Performing the operation x(t) = <X(t)> +
    sum(sqrt(lambda_i)*phi_i(t)*Y_i)
%}
x = dataKLE.meanX(t);


nY = size(Y,2);
if nY > 1
    x = repmat(x, 1, nY);
end

% Matrix goes fast?
x = x + phiMatrix*(dataKLE.lambda.^(0.5) .* Y);
end