function xTotDot = DynPlusJacobian(t, xTot, fController, const, theta)
    %{
        Calculates the jacobian of state x with respect to Cm_alpha.
    %}
    x = xTot(1:5);
    Sk = xTot(6:end);
        
    % Polynomial CM_alpha
    Cm_Alpha = 0;
    nTheta = numel(theta);
    for i = 1:nTheta
        Cm_Alpha = Cm_Alpha + theta(i)*t.^(i-1);
    end

    % Vanilla Dynamics
    [u,~] = fController(t, x, const);
    xDot = getDynamics(x, u, const, t);

    % Jacobian integration
    Jinv = 1/const.gravimetrics.Iyy;
    refArea = const.aero.refArea;
    cBar = const.aero.refLength;
    alpha = calcAlpha(x);
    qBar = 1/2 * getAtmo(x(5), const) * x(1)^2;

    B = [0; 0; Jinv*qBar*refArea*cBar*alpha; 0; 0]; %df_dCMalpha
    A = AmatrixWrapper(x, const, Cm_Alpha, u);
    
    Sdot = A*Sk + B;

    xTotDot = [xDot; Sdot];

end