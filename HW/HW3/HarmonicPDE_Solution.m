function [u] = HarmonicPDE_Solution(uBC, K, h)
%{
    Solves the 1D heat diffusion PDE for HW3 using flux form:
    d/dx(K(x) * d/dx(u(x)) = 1

    Uses harmonic average of K.

    Inputs:
        uBC = 2x1 vector of boundary conditions at x = 0, x = N
        K = N x 1 vector of thermal conductivity values
        h = finite discretization step size 
%}

    N = size(K,1);
    u0 = uBC(1);
    uN = uBC(2);


    % Build A matrix
    %% TODO - might need to truncate this matrix to avoid issues?
    Ki = K(1:N-1);
    Kip1 = K(2:N);
    Kit1 = K(1:N-1); % this feels wrong?
    K_harmonicMeanPlus = 2 * Ki .* Kip1 ./ (Ki + Kip1);
    K_harmonicMeanMinus = 2 * Ki .* Kit1 ./ (Ki + Kit1);
    diagTerms = 
    A = zeros(N,N);


    % Build b vector
    b = [u0; ones(N-2,1)*h^2; uN];

    % Solve System
    u = A\b; 

    % Force BC to avoid any numerical precision issues
    u(0) = uBC(1);
    u(N) = uN;
end