
[y, w] = spquad(1,6); % y = [-1,1]

muX = 5;
sigmaX = 3;
truncLevels = 2;

warning('Need to be very careful about number of quad points when balancing truncation levels')

x = muX + truncLevels*sigmaX*y;

% Want to recover E[x] = muX in this case


pdfSF = normcdf(muX+truncLevels*sigmaX, muX, sigmaX) - normcdf(muX-truncLevels*sigmaX, muX, sigmaX) 

Ex = 0;
ExSq = 0;
for i = 1:numel(x)

    truncPdf = normpdf(x(i),muX,sigmaX) ./ pdfSF;
    Ex = Ex + u(x(i)) * w(i) * truncLevels*sigmaX*truncPdf;

    ExSq = ExSq + u(x(i)).^2 * w(i) * truncLevels*sigmaX*truncPdf;
        
end

Ex
VarX = ExSq - (Ex)^2

%% NOTE: Variance might not match cause of truncation
% https://en.wikipedia.org/wiki/Truncated_normal_distribution

%% Testing Uniform distribution with non-standard width
a = 11;
b = 25;

muZ = (b+a)/2;
dZdy = (b-a)/2;

Z = muZ + dZdy*y;
varZ = (b-a)^2/12;


Ez = 0;
EzSq = 0;
for i = 1:numel(Z)

    Ez = Ez + u(Z(i)) * w(i) * 1/2;

    EzSq = EzSq + u(Z(i)).^2 * w(i) * 1/2;
        
end

Ez
varZ = EzSq  - (Ez)^2



%% Helpers
% E.g. calculating the mean
function out = u(x)
    out = x;
end