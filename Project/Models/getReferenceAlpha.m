function [alpha, q] = getReferenceAlpha(t)
    %{
        Generates a sine sweep in alpha with amplitude of 10 degrees and
        frequency of 0.1Hz.
    %}  
    f = 1/10;
    omega = 2*pi*f;
    A = deg2rad(10);
    alpha = A * sin(omega*t);
    q = -A*omega*cos(omega*t);    
end