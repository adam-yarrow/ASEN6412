function [dataKLE] = generateKLE_Struct(lambda, omega, meanX, a, d)
%{
    Builds a structure to hold the reusuable componets of a KLE
    approximation of a random variable.
%}
    
    dataKLE = struct();
    dataKLE.lambda = lambda;
    dataKLE.omega = omega;
    dataKLE.meanX = meanX;
    dataKLE.a = a;
    dataKLE.d = d;
    dataKLE.phiVector = cell(d, 1);

    for i = 1:d
        dataKLE.phiVector{i} = evalEigenfunctionKLE(omega(i), a, i);
    end
    
end