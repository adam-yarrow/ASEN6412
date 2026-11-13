function h = CreateMeasEstimates(theta, const, controller)
%{
    Runs model for a given set of theta
%}
    const.aero.polyCM_alpha.coeff = theta;
    [simData] = Simulation(controller, const);

    h = pagemtimes(const.measModel.C, simData.x);
    h = reshape(h,[],1);    
end