function [muU, varU] = StochasticColloc(Y_Tensors, W_Tensor, Nx, uBC, ...
    dataKLE, dx, integrationSF)
    %{
        Y_Tensors = cell array (NumLevels x 1) of KLE random variables for
        each stochastic collocation level.
    
        W_Tensor = cell array (NumLevels x 1) of weights for stochastic
        collocation.

        Nx = number of discretization points

        uBC = 2x1 boundary conditions

        dataKLE = KLE data structure

        dx = discretization grid size

        integrationSF = scale factor for converting integration into
        expectation operator.

        RETURNS:
            muU = mean of u(x) 
            varU = variance of u(x)
    %}

% TODO - not super efficient, could replace with trace(W^T*U(y1,y2))
% OR Reuse function evals
% Evaluate Expected Values


    NumLevels = numel(W_Tensor);
    
    muU = NaN(Nx, NumLevels);
    varU = NaN(Nx, NumLevels);
    
    KipHalf_vs_x = NaN(Nx,1); % padded with NaNs
    KitHalf_vs_x = NaN(Nx,1);
    
    for iLevel = 1:NumLevels   
        NquadPoints = size(Y_Tensors{iLevel},1);
        Iquad_Mu = 0; % E[u]
        Iquad_uSq = 0; % E[u^2]
    
        Y_CC = Y_Tensors{iLevel};
        W_CC = W_Tensor{iLevel};   
    
        for i = 1:NquadPoints        
            KitHalf_vs_x(2:end) = getRV_KLE_Fast(dataKLE.x_itHalf, dataKLE, ...
                                        dataKLE.phiMatrix, Y_CC(i,:)'); 
            KipHalf_vs_x(1:end-1) = KitHalf_vs_x(2:end); % Symmetry condition - saves compute
            
            uY = PDE_Solution(uBC, KipHalf_vs_x, KitHalf_vs_x, dx);
       
            Iquad_Mu = Iquad_Mu + uY * W_CC(i);
            Iquad_uSq = Iquad_uSq + uY .^2 * W_CC(i);
        end
    
        UcrossMoment = Iquad_uSq * integrationSF;
        muU(:, iLevel) = Iquad_Mu * integrationSF;
        varU(:, iLevel) = UcrossMoment - muU(:, iLevel).^2;
    end


end