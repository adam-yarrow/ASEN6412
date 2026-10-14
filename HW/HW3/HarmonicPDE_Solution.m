function [u] = HarmonicPDE_Solution(uBC, KipHalf, KitHalf, h)
%{
    Solves the 1D heat diffusion PDE for HW3 using flux form:
    d/dx(K(x) * d/dx(u(x)) = 1

    Inputs:
        uBC = 2x1 vector of boundary conditions at x = 0, x = N
        KipHalf = N x 1 vector of thermal conductivity values evaluated at
             i+1/2 (valid for i = 1:N-1)
        KitHalf = N x 1 vector of thermal conductivity values evaluated at
             i-1/2 (valid for i = 2:N)
        h = finite discretization step size 
%}

    N = length(KipHalf);
    u0 = uBC(1);
    uN = uBC(2);

    % Build A Matrix (sparsely)
    Aii = [1; KipHalf(2:N-1) + KitHalf(2:N-1); 1]; % Range goes from i = 2 to i = N-1
    lowerDiag = [-KitHalf(2:N-1); 0; 0];
    upperDiag = [0; 0; -KipHalf(2:N-1)];
    A = spdiags([lowerDiag, Aii, upperDiag], -1:1, N, N);
    
    % Build b vector
    b = [u0; ones(N-2,1)*h^2; uN];

    % Solve System
    u = A\b; 

    % Force BC to avoid any numerical precision issues
    u(1) = uBC(1);
    u(N) = uN;
end