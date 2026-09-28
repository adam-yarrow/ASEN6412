function [x] = getRV_KLE(t, dataKLE, Y)
 %{
    Returns the random value of the stochastic process at
    location/time t.

    Performing the operation x(t) = <X(t)> +
    sum(sqrt(lambda_i)*phi_i(t)*Y_i)
%}
x = dataKLE.meanX(t);

nY = size(Y,2);
if nY >1
    x = repmat(x, 1, nY);
end

% Translate t into domain of eigen function
tTilde = t - dataKLE.a;

% TODO - try to make this into a vector operation??? e.g. inner
% product
for i = 1:dataKLE.d
    x = x + sqrt(dataKLE.lambda(i)).*dataKLE.phiVector{i}(tTilde)*Y(i,:);
end
end