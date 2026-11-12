%{
    Need to test stochastic collocation
%}

%% u(x) = x^2, and x ~ N(mu, sigma) with caps at +/- K-sigma about the mean
mu = 2.12;
sigma = 2.13/10;
K = 3;

%% Analytical Solution - equivalent to E[u] = E[X^2]
% From wikipedia
beta = K;
alpha = -K;
Z = normcdf(beta) - normcdf(alpha);

varX = sigma^2*(1 - (beta*normpdf(beta) - alpha*normpdf(alpha))/Z ...
            - (normpdf(alpha) - normpdf(beta))/Z);
meanX = mu + (normpdf(alpha) - normpdf(beta))/Z*sigma;

Eu_Analytical = varX + meanX^2;

%% SC of solution
d = 1;
q = 2;
[Z, W] = spquad(d, q);

Y = mu + K*sigma*Z;

Wcorrected = K * sigma * W .*truncNormalPdf(Y, mu, sigma, K);
Wsum = sum(Wcorrected);

Eu = 0;
for i = 1:numel(W)
    uVal = u(Y(i));

    Eu = Eu + uVal*Wcorrected(i);    
end

fprintf('q = %i\n', q);
fprintf('Eu Approx = %3.3f\n', Eu);
fprintf('Eu Exact = %3.3f\n', Eu_Analytical)
fprintf('Eu Error (percent) = %3.3f\n', (Eu - Eu_Analytical)/Eu_Analytical * 100)
fprintf('sum(W) = %1.3f\n', Wsum);


%% Claude suggestion

[Zq, ~] = spquad(1, q);
zn = sort(Zq(:));                          % nodes in [-1,1]
n  = numel(zn);
pz = @(z) K*sigma .* truncNormalPdf(mu + K*sigma*z, mu, sigma, K);  % density in z
m  = arrayfun(@(k) integral(@(z) z.^k .* pz(z), -1, 1, ...
        'AbsTol',1e-14,'RelTol',1e-13), 0:n-1).';                   % moments
Wex = (zn(:).^(0:n-1)).' \ m;              % solves sum_i w_i z_i^k = m_k
Y   = mu + K*sigma*zn;
EuClaude  = sum(u(Y) .* Wex);


%% Try bivariate distribution after this


function out = u(x)
    out = x.^2;
end


