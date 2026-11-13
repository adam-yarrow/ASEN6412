function A = AmatrixWrapper(xk, const, CM_alpha, deltaE)
    %{
        Wraps up symbolic math jacobian to linearize dynamics. Allows for
        Cm_Alpha to be provided by the user (to account for varying
        CM_alpha).

        NOTE: Jacobian assumes density is invariant w.r.t. height
        (reasonable for small height changes).
    %}
    rho = getAtmo(xk(5), const);
    cAero = const.aero;
    cGrav = const.gravimetrics;
    A = compute_A_matrix(xk,rho,cGrav.m,cGrav.Iyy,const.g,...
        const.aero.refLength,const.aero.refArea,...
        deltaE,...
        cAero.CL_alpha,cAero.CL_deltaE,cAero.CL_0,...
        cAero.CD_alphaSq,cAero.CD_alpha,cAero.CD_deltaE,cAero.CD_deltaESq,cAero.CD_0,...
        cAero.CM_alphaSq,CM_alpha,cAero.CM_deltaE,cAero.CM_0);

end