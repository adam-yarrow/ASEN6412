%% Generating Jacobian of nonlinear dynamics symbolically
% NOTE: assuming density is constant and not a function of altitude

syms V gamma q theta h
x = [V; gamma; q; theta; h];

syms rho m J g cBar S deltaE 
syms CL_alpha CL_deltaE CL_0
syms CD_alphaSq CD_alpha CD_deltaE CD_deltaESq CD_0
syms CM_alphaSq CM_alpha CM_deltaE CM_0

alpha = theta - gamma;
qBar = 1/2*rho*V^2;

CL = CL_alpha*alpha + CL_deltaE*deltaE + CL_0;

CD = CD_alphaSq*alpha.^2 + CD_alpha*alpha + ...
    CD_deltaE*deltaE + CD_deltaESq*deltaE.^2 + ...
    CD_0;

CM = CM_alphaSq*alpha.^2 + CM_alpha*alpha + CM_deltaE*deltaE + ...
     CM_0;

D = qBar.*S.*CD;
L = qBar.*S.*CL;
M = qBar.*S.*cBar.*CM;
    

xDot(1) = (1/m)*(-D) - g*sin(gamma);   
xDot(2) = (1./(m*V)).*(L) - g./V .* cos(gamma);
xDot(3) = M/J;
xDot(4) = q;
xDot(5) = V.*sin(gamma); 

dFdx = jacobian(xDot, x);

disp('Symbolic Jacobian:');
disp(dFdx);

% Generate the callable numerical function
matlabFunction(dFdx, 'File', 'compute_A_matrix', ...
    'Vars', {x, rho m J g cBar S deltaE CL_alpha CL_deltaE CL_0 CD_alphaSq CD_alpha CD_deltaE CD_deltaESq CD_0 CM_alphaSq CM_alpha CM_deltaE CM_0});
% Inputs to the function : debris_state_k, gps_state_k, rx_state_k, fT_val, c_val