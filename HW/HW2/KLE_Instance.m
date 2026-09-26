classdef KLE_Instance < handle
    %{
        This class defines an instance of a stochastic process X(omega,t) that is
        approximated by a truncated KL Expansion.
    
        This allows the user to call the class multiple times to get the
        value of the approximate process at different spatial/temporal 
        locations for a fixed random draw.

        TODO - this class might have prohibitive overhead when dealing with
        100,000 plus instantiations, in which case a class passing in the
        data struct might be nicer. TBD see how it goes
    %}  

    properties
        Y; % Random variable draws for this instance (d x 1)
        lambda; % Eigenvalues for KLE expansion (n x 1) NOTE: n >= d
        omega; % Frequencies for KLE expansion (n x 1) 
        meanX; % Mean of the stochastic process that has been expanded as a 
                % function of time or space
        a; % half domain of eigenfunction (e.g. -a < t_tilde < a)
        d; % Truncation order
        phiVector; % A vector of phi functions scaled by 
                            % sqrt(lambda) (dx1)

    end

    methods
        function obj = KLE_Instance(Y, lambda, omega, meanX, a)
            % Initialize variables
            obj.Y = Y;
            obj.lambda = lambda;
            obj.omega = omega;
            obj.meanX = meanX;
            obj.a = a;

            obj.d = length(Y);

            obj.phiVector = cell(obj.d,1);
            for i = 1:obj.d
                obj.phiVector{i} = evalEigenfunctionKLE(obj.omega(i), ...
                                        obj.a, i);
            end
        end

        function x = getRV(obj,t)
            %{
                Returns the random value of the stochastic process at
                location/time t.

                Performing the operation x(t) = <X(t)> +
                sum(sqrt(lambda_i)*phi_i(t)*Y_i)
            %}
            x = obj.meanX(t);

            % Translate t into domain of eigen function
            tTilde = t - obj.a;

            % TODO - try to make this into a vector operation??? e.g. inner
            % product
            for i = 1:obj.d
                x = x + sqrt(obj.lambda(i)).*obj.phiVector{i}(tTilde).*obj.Y(i);
            end

        end
    end
end