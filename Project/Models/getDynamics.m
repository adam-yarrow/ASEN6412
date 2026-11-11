function [xDot] = getDynamics(x, u, const, t)
%{
    Gets xDot = f(x,u)
    x = [V, gamma, q, theta, h]
%}   
    deltaE = u(1);

    % States and Derived params
    V = x(1);
    gamma = x(2);
    q = x(3);
    theta = x(4);
    h = x(5);
    alpha = calcAlpha(x);

    % Sys Params
    m = const.gravimetrics.m;
    J = const.gravimetrics.Iyy;
    g = const.g;

    cBar = const.aero.refLength;
    S = const.aero.refArea;

    % Get Dynamic Pressure
    [rho, ~] = getAtmo(h, const);
    qBar = 0.5.* rho .* V.^2; % Pa

    % Get External Forces
    [CL, CD, CM, CMq] = getAeroCoeffs(alpha, deltaE, const, t, x);

    D = qBar.*S.*CD;
    L = qBar.*S.*CL;
    M = qBar.*S.*cBar.*(CM + (cBar./(2*V)).*CMq);
    
    % Dynamics
    xDot = zeros(const.nStates,1);
    xDot(1) = (1/m)*(-D) - g*sin(gamma);   
    xDot(2) = (1./(m*V)).*(L) - g./V .* cos(gamma);
    xDot(3) = M/J;
    xDot(4) = q;
    xDot(5) = V.*sin(gamma);    
end