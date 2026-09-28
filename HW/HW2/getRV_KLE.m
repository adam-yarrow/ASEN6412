function [x] = getRV_KLE(t, dataKLE, Y)
 %{
    Returns the random value of the stochastic process at
    location/time t.

    Performing the operation x(t) = <X(t)> +
    sum(sqrt(lambda_i)*phi_i(t)*Y_i)
%}
x = dataKLE.meanX(t);

% Translate t into domain of eigen function
tTilde = t - dataKLE.a;

% TODO - try to make this into a vector operation??? e.g. inner
% product
for i = 1:obj.d
    x = x + sqrt(dataKLE.lambda(i)).*dataKLE.phiVector{i}(tTilde).*Y(i);
end
end