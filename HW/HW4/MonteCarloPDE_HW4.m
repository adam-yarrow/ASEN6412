function [u, uStats] = MonteCarloPDE_HW4(x, uBC, Nmc, dataKLE, Y)
    %{
    Function to evaluation the MC approximation of the PDE for a fixed grid
    size.

    x = grid to evaluate over
    Nmc = Number of monte-carlo runs
    uBC = 1 x 2 vector containing the boundary conditions for each run
    dataKLE = KLE data structure
    Y = random variables used for KLE expansion of K(x,w), size d x
    Nmc
    fPlot = Flag to enable/disable plotting
    %}

    Nx = numel(x);

    %% Solve MC
    u = NaN(Nx,Nmc);
    KipHalf_vs_x = NaN(Nx,1); % padded with NaNs
    KitHalf_vs_x = NaN(Nx,1);

    dx = x(2) - x(1);

    % Run MC
    for i = 1:Nmc
        % Sampling K(x) - doing inside loop to avoid overflowing memory for
        % large Nx and Nmc - at cost of more compute time
        KitHalf_vs_x(2:end) = getRV_KLE_Fast(dataKLE.x_itHalf, dataKLE, dataKLE.phiMatrix, Y(:,i)); 
        KipHalf_vs_x(1:end-1) = KitHalf_vs_x(2:end); % Symmetry condition - saves compute
        
        % Solving
        u(:,i) = PDE_Solution(uBC, KipHalf_vs_x, KitHalf_vs_x, dx);
    end

    %% Generate Statistics
    xHalfIdx = find(x == 0.5); %% Brittle - fix
    uStats = struct();
    uStats.mu = mean(u, 2);
    uStats.var = var(u, [], 2);
    uStats.muAtXhalf = uStats.mu(xHalfIdx, :);
    uStats.varAtXhalf = uStats.var(xHalfIdx, :);

end
