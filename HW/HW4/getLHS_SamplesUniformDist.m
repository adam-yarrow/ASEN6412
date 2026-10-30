function Zi = getLHS_SamplesUniformDist(N,a,b)
    %{
        Generates LHS samples assuming variables are iid uniform
        distributions.
        
        N = number of samples.
        a = d x 1 lower bounds for uniform distribution
        b = d x 1 upper bounds for uniform distribution
    %}

    d = numel(a);

    Zi = NaN(d, N);

    for i = 1:d
        % Splitting uniform domain into N equal probability domains Zi
        divisions = linspace(a(i), b(i), N+1);
        Zlower = divisions(1:end-1);
        Zupper = divisions(2:end);

        % Draw 1 sample from each Z division
        Zi(i,:) = Zlower + (Zupper - Zlower) .* rand(1,N);

        % Randomize ordering for combining dimensions
        Zi(i,:) = Zi(i, randperm(N));
    end
end