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

    % force dimension of K to be column vector
    if size(K,1) < size(K,2)
        K = K';
    end

    N = length(K);
    u0 = uBC(1);
    uN = uBC(2);

%% TODO - pass in Kiphalf and KitHalf from K(x,w) solver....
    warning('TO FIX')
    % Build A diagonals
    KipHalf = (2 * K(1:N-1) .* K(2:N)) ./ (K(1:N-1) + K(2:N)); % Domain: 1 <= i <= N-1
    KitHalf = (2 * K(2:N) .* K(1:N-1)) ./ (K(2:N) + K(1:N-1)); % Domain: 2 <= i <= N; This is identical to KipHalf - can i just remove?

    Aii = [1; KipHalf(2:end) + KitHalf(1:end-1); 1];
    A = diag(Aii);

    % Add off diagonals
    A = A + ...
        diag([-KitHalf(1:N-2); 0], -1) + ...
        diag([0; -KipHalf(2:N-1)], +1);

    % Build b vector
    b = [u0; ones(N-2,1)*h^2; uN];

    % Solve System
    u = A\b; 

    % Force BC to avoid any numerical precision issues
    u(1) = uBC(1);
    u(N) = uN;
end