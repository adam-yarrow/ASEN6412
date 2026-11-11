%{
    Need to test stochastic collocation
%}

%% u(x) = x^2, and x ~ N(mu, sigma) with caps at +/- K-sigma about the mean
mu = 2.12;
sigma = 2.13/10;
K = 3;

% Analytical Solution
L = 2*normcdf(K) - 1;
Eu_analytical = 1/L*((mu^2+sigma^2)*(2*normcdf(K) - 1) - 2*sigma^2*K*normpdf(K));


%% Use analytical formulat from wikipedia too.
% Symmetric bounds so can solve this


%% TODO SC