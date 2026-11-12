function [x] = getRV_KLE(t, dataKLE, Y)
 %{
    Returns the random value of the stochastic process at
    location/time t.

    Performing the operation x(t) = <X(t)> +
    sum(sqrt(lambda_i)*phi_i(t)*Y_i)
%}
x = dataKLE.meanX(t);
nT = numel(t);

nY = size(Y,2);
if nY > 1
    x = repmat(x, 1, nY);
end

% Translate t into domain of eigen function
tTilde = t - dataKLE.a;

% Pre-make eigenfunctions
phiMatrix = NaN(nT, nY);
for i = 1:dataKLE.d
    phiMatrix(:,i) = evalEigenfunctionKLE(dataKLE.omega(i), dataKLE.a, i, tTilde);
end
lambda = dataKLE.lambda(:);

% Matrix goes fast?
x = x + phiMatrix*(lambda.^(0.5) .* Y);
end