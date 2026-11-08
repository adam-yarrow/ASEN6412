function [fx] = truncNormalPdf(x, mu, sigma, truncSigmaLevel)
    %{
        Calculates the truncated gaussian density for a specified value.

        Reference:
        https://en.wikipedia.org/wiki/Truncated_normal_distribution
    %}
    
    % SF = FX(upper bound) - FX(lower bound)
    pdfSF = normcdf(mu + truncSigmaLevel*sigma, mu, sigma) - ...
            normcdf(mu - truncSigmaLevel*sigma, mu, sigma);

    fx = normpdf(x, mu, sigma) / pdfSF;
end