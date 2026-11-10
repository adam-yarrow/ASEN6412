%% Frequency Domain Fuck around
%{
    Look at MC response of q, theta, gamma and alpha in frequency domain
    from T = 3 seconds onwards. A, phase, frequency mean and variance of interest
%}

validIdx = find(times >= 3);

stateIdx = 3; % q to test

q = squeeze(rawData_MC.x(stateIdx,validIdx,:)) * rad2deg(1);

% detrend
qNoTrend = detrend(q,"linear"); % works on each column independently (e.g. each MC series)

qTrend = q - qNoTrend;
qTrendMean = mean(qTrend,2); % Mean of all trends at each time

% Convert to frequency domain
[f, A, Xhalf, phi] = onesidedFFT(qNoTrend(:,1), 1/const.dT, true, 1E-6);

%% NOTE: Issue i have is my signal decays in amplitude over time, so FFT will lose some of this info


%% TODO - fix up this function to be properly useful
%% TODO - consider a spectrogram?
function [f, A, Xhalf, phi] = onesidedFFT(x, fs, fPlot, tol)
    %{
        One-sided amplitude and phase spectrum of a real time series.

        Inputs:
            x      - real vector (row or column)
            fs     - sample rate in Hz
            doPlot - true to plot amplitude and phase
            tol    - phase is hidden (NaN) for bins whose amplitude is below
                     tol*max(A), since phase is meaningless there 

        Outputs:
            f     - frequency vector, 0 to Nyquist (Hz), column
            A     - one-sided amplitude spectrum (same units as x)
            Xhalf - one-sided complex spectrum, scaled by 1/N
            phi   - phase in degrees, relative to a cosine, NaN where masked
    %}

    x = x(:);
    N = numel(x);

    X = fft(x) / N;
    nHalf = floor(N/2) + 1;
    Xhalf = X(1:nHalf);

    if mod(N,2) == 0
        Xhalf(2:end-1) = 2*Xhalf(2:end-1);   % even N: DC and Nyquist unpaired
    else
        Xhalf(2:end)   = 2*Xhalf(2:end);     % odd N: only DC unpaired
    end

    A   = abs(Xhalf);
    phi = rad2deg(angle(Xhalf));
    phi(A < tol*max(A)) = NaN;               % hide noise-level phase

    f = (0:nHalf-1).' * fs / N;

    if fPlot
        figure('Name','One-sided spectrum');
        tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        ax1 = nexttile(tl);
        plot(ax1, f, A, 'LineWidth', 1.2);
        grid(ax1,'on'); ylabel(ax1,'Amplitude');

        ax2 = nexttile(tl);
        plot(ax2, f, phi, '.-');
        grid(ax2,'on'); ylabel(ax2,'Phase (deg)'); xlabel(ax2,'Frequency (Hz)');
        ylim(ax2, [-180 180]); yticks(ax2, -180:90:180);

        linkaxes([ax1 ax2], 'x');
    end
end