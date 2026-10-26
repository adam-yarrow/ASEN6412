function [alpha, q] = getReferenceAlpha(t)
    f = 1/10;
    omega = 2*pi*f;
    A = deg2rad(10);
    alpha = A * sin(omega*t);
    q = -A*omega*cos(omega*t);    
end