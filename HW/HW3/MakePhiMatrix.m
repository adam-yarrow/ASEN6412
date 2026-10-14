function [phiMatrix] = MakePhiMatrix(t, dataKLE)
nT = numel(t);

% Translate t into domain of eigen function
tTilde = t - dataKLE.a;

phiMatrix = NaN(nT, dataKLE.d);
for i = 1:dataKLE.d
    phiMatrix(:,i) = evalEigenfunctionKLE(dataKLE.omega(i), dataKLE.a, i, tTilde);
end
end