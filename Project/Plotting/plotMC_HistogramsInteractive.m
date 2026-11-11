function plotMC_HistogramsInteractive(const, rawData, nSamplesPerBin, initialTime)
    %{
        Interactive histograms of the MC state response.
        Slider scrubs through time; the edit box sets the time step
        between frames (in the same units as rawData.simData{1}.times).

        AI GENERATED (there's no way im that good at making plots)
    %}
    times  = rawData.simData{1}.times(:);
    dtMin  = median(diff(times));            % native sample spacing
    dt     = dtMin;
    if nargin < 4, initialTime = times(1); end
    if nargin < 3, nSamplesPerBin = 50; end

    nPlots = const.nStates + 1;              % states + alpha
    frameTimes = []; frameIdx = [];          % set in rebuildFrames

    % ---------- figure and axes ----------
    fig = figure('Name','MC State Variable Histograms', ...
                 'Units','normalized','Position',[0.05 0.08 0.9 0.82]);
    tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0 0.10 1 0.90];      % leaves a strip at the bottom for controls


    h = gobjects(nPlots,1);
    for k = 1:nPlots
        ax = nexttile(tl);
        h(k) = histogram(ax, nan, 'Normalization','probability');
        grid(ax,'on');
        ylabel(ax,'Probability');
        if k <= const.nStates
            xlabel(ax, sprintf('%s (%s)', const.stateNames{k}, ...
                               const.statePlottingUnits{k}));
        else
            xlabel(ax, 'Alpha (deg)');
        end
    end

    % ---------- controls ----------
    playHz = 1;                                   % playback rate, frames per second

    btnPlay = uicontrol(fig,'Style','pushbutton','Units','normalized', ...
        'Position',[0.02 0.03 0.06 0.04],'String','Play', ...
        'Callback',@onPlayPressed);
    uicontrol(fig,'Style','text','Units','normalized', ...
        'Position',[0.09 0.03 0.04 0.04],'String','Time:', ...
        'HorizontalAlignment','right');
    sld = uicontrol(fig,'Style','slider','Units','normalized', ...
        'Position',[0.14 0.03 0.50 0.04],'Min',1,'Max',2,'Value',1);
    lblTime = uicontrol(fig,'Style','text','Units','normalized', ...
        'Position',[0.65 0.03 0.11 0.04],'String','', ...
        'HorizontalAlignment','left');
    uicontrol(fig,'Style','text','Units','normalized', ...
        'Position',[0.77 0.03 0.09 0.04],'String','Time step:', ...
        'HorizontalAlignment','right');
    edtDt = uicontrol(fig,'Style','edit','Units','normalized', ...
        'Position',[0.87 0.03 0.08 0.04],'String',num2str(dt), ...
        'Callback',@onStepChanged);

    tmr = timer('ExecutionMode','fixedRate','BusyMode','drop', ...
                'Period',1/playHz,'TimerFcn',@(~,~) onTick());

    addlistener(sld,'Value','PostSet',@(~,~) onSliderMoved());  % live scrubbing
    fig.WindowKeyPressFcn = @onKey;                             % left/right arrows
    fig.CloseRequestFcn   = @onClose;                           % clean up timer

    rebuildFrames(initialTime, nSamplesPerBin);                             

    % ================= nested functions =================
    function rebuildFrames(tKeep, nSamplesPerBin)
        frameTimes = (times(1):dt:times(end))';
        if numel(frameTimes) < 2, frameTimes = times([1 end]); end
        frameIdx = interp1(times, 1:numel(times), frameTimes, 'nearest');
        n = numel(frameTimes);
        [~,k0] = min(abs(frameTimes - tKeep));
        sld.Min = 1; sld.Max = n; sld.Value = k0;
        sld.SliderStep = [1/(n-1), min(1,10/(n-1))];
        drawFrame(k0, nSamplesPerBin);
    end

    function onSliderMoved()
        drawFrame(round(sld.Value), nSamplesPerBin);
    end

    function onStepChanged(~,~)
        newDt = str2double(edtDt.String);
        if isnan(newDt) || newDt <= 0
            edtDt.String = num2str(dt);  return
        end
        dt = max(newDt, dtMin);          % can't be finer than the data
        edtDt.String = num2str(dt);
        rebuildFrames(frameTimes(round(sld.Value)), nSamplesPerBin);
    end

    function onKey(~,evt)
        switch evt.Key
            case 'rightarrow', sld.Value = min(sld.Max, round(sld.Value)+1);
            case 'leftarrow',  sld.Value = max(sld.Min, round(sld.Value)-1);
        end
    end

    function drawFrame(k, nSamplesPerBin)
        ti = frameIdx(k);
        for i = 1:const.nStates
            data = const.statePlottingSF(i) * squeeze(rawData.x(i,ti,:));
            setHist(h(i), data, nSamplesPerBin);
        end
        setHist(h(nPlots), rad2deg(1) * squeeze(rawData.alpha(ti,:))', nSamplesPerBin);
        lblTime.String = sprintf('T = %.2f', times(ti));
        title(tl, sprintf('MC State Variable Histograms at T = %.2f', times(ti)));
        drawnow limitrate
    end

    function setHist(hh, data, nSamplesPerBin)
        hh.Data    = data(:);
        hh.NumBins = max(1, round(numel(data)/nSamplesPerBin));
    end

    function onPlayPressed(~,~)
        if strcmp(tmr.Running,'on')
            stopPlayback();
        else
            if round(sld.Value) >= sld.Max, sld.Value = sld.Min; end  % restart from the beginning
            btnPlay.String = 'Pause';
            start(tmr);
        end
    end

    function onTick()
        if ~isvalid(fig), return, end
        k = round(sld.Value);
        if k >= sld.Max
            stopPlayback();                 % reached the last frame
        else
            sld.Value = k + 1;              % listener triggers drawFrame
        end
    end

    function stopPlayback()
        if isvalid(tmr) && strcmp(tmr.Running,'on'), stop(tmr); end
        btnPlay.String = 'Play';
    end

    function onClose(~,~)
        if isvalid(tmr), stop(tmr); delete(tmr); end
        delete(fig);
    end
end