function fn_out = ephysLoosePatchSpikeFinder(fpath, varargin)
% ephysLoosePatchSpikeFinder finds spikes in loose patch data
% fn_out = ephysLoosePatchSpikeFinder(fpath, varargin)

%% Parse
if nargin < 2
    varargin = {};
    if nargin < 1
        fpath = '';
    end
end

% Debug
% fpath = 'E:\ephys\stephen\May 2021 loose AVPV TH Stephen\210513a loose AVPV TH Stephen\210513a loose AVPV TH Stephen_0001.abf';

p = inputParser;

% Data handling parameters
addOptional(p, 'defaultpath', '\\anastasia\data\ephys\stephen\*.abf');
addOptional(p, 'defaultpathmat', '\\anastasia\data\ephys\stephen\*_preprocess.mat');
addOptional(p, 'useprevparams', true); % Use previous parameters
addOptional(p, 'loadmat', false); % Load matlab mat files
addOptional(p, 'segment', []); % Analyzing different parts of the trace

% Analysis parameters
addOptional(p, 'channel', 1);
addOptional(p, 'invert', true);
addOptional(p, 'subtractmovingmedian', true); % Subtract moving median or statis median
addOptional(p, 'movingmedianwin', 200) % Number of points in the moving median window
addOptional(p, 'defaultprom', 10); % Default prominence threshold
addOptional(p, 'defaultcurrent', 5); % Default current threshold
addOptional(p, 'defaultwidth', 5); % Default wid threshold
addOptional(p, 'defaultdis', 2); % default distance in ms
addOptional(p, 'usedatawindow', true); % Use a window of data

% Display parameters
addOptional(p, 'pos', [50 200 1800 500]); % Figure position
addOptional(p, 'spikewindow', [-2 2]); % Window to display spikes, in ms
addOptional(p, 'nbins', 20);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% IO
% Load
if isempty(fpath) && p.loadmat
    [fn, fp] = uigetfile(p.defaultpathmat);
    [~, fn, ~] = fileparts(fn);
    fn = fn(1:strfind(fn, '_preprocess')-1);
    
    % Full filename
    fpath = fullfile(fp, sprintf('%s%s', fn, '.abf'));
elseif isempty(fpath)
    [fn, fp] = uigetfile(p.defaultpath);
    [~, fn, ext] = fileparts(fn);
    
    % Full filename
    fpath = fullfile(fp, sprintf('%s%s', fn, ext));
else
    [fp, fn, ~] = fileparts(fpath);
end

% Full outputname
if isempty(p.segment)
    fn_out = fullfile(fp, sprintf('%s_preprocess.mat', fn));
else
    fn_out = fullfile(fp, sprintf('%s_preprocess_s%i.mat', fn, p.segment));
end

% Load
[d, si, dataparams] = abfload(fpath);

% Read parameters
if exist(fn_out, 'file') && p.useprevparams
    params = load(fn_out, '-mat');
    newparams = false;
    disp('Using previous parameters')
    
    if isfield(params, 'ws')
        reload = false;
        ws = params.ws;
    else
        reload = true;
    end
else
    newparams = true;
end

%% Data window (to remove time periods when a cell is lost)
if p.usedatawindow && newparams
    hdatawin = figure('Position', p.pos, 'Name', 'Choose data window');
    plot(d(:,p.channel));
    
    % Choose region if needed
    boxui = imrect(gca);
    userbox = wait(boxui);
    delete(boxui)
    datawin = round([userbox(1), userbox(1) + userbox(3)]);
    datawin(1) = max(1, datawin(1));
    datawin(2) = min(length(d), datawin(2));

    % Calculate data that are in the box
    d = d(datawin(1) : datawin(2))';
    close(hdatawin)
elseif p.usedatawindow && ~newparams
    datawin = params.datawin;
    d = d(datawin(1) : datawin(2))';
else
    datawin = [1 length(d)];
end

%% Preprocess
% Get channel
d = d(:, p.channel);
l = length(d);

% Sampling rate and times
fs = 1000/si; % in ms-1
tmax = l/fs/1000;
t0 = 1/fs/1000;
tvec = t0 : 1/fs/1000 : tmax;

% distance
dis = p.defaultdis * fs;

% Apply inversions and moving medians
if p.invert 
    d2 = -d;
else
    d2 = d;
end
if p.subtractmovingmedian
    movmed = movmedian(d2, p.movingmedianwin);
    d2 = d2 - movmed;
else
    med = median(d2);
    d2 = d2 - med;
end

% Thresholds
if newparams
    promt = p.defaultprom;
    currt = p.defaultcurrent;
    widtht = p.defaultwidth;
else
    promt = params.promt;
    currt = params.currt;
    if isfield(params, 'widtht')
        widtht = params.widtht;
    else
        widtht = p.defaultwidth;
    end
end

% Inds for spikes
inds = p.spikewindow(1)*fs : p.spikewindow(2)*fs;

%% Figure
hfig = figure('Position', p.pos);

% Spike window
subplot(1,6,1:4)
if newparams
    [~, locs, ws, proms] = findpeaks(d2, 'MinPeakProminence', promt,...
        'MinPeakHeight', currt, 'MinPeakDistance', dis, 'MinPeakWidth', widtht);
elseif reload
    [~, locs, ws, proms] = findpeaks(d2, 'MinPeakProminence', promt,...
        'MinPeakHeight', currt, 'MinPeakDistance', dis, 'MinPeakWidth', widtht);
else
    locs = params.locs;
    proms = params.proms;
end
hold on
plot(tvec, d(:,1))
hloc = plot(tvec(locs), d(locs), 'r.');
if p.subtractmovingmedian
    if p.invert
        hcurr = plot(tvec, -movmed - currt, 'Color', [0.8 0.8 0.8]);
    else
        hcurr = plot(tvec, movmed + currt, 'Color', [0.8 0.8 0.8]);
    end
else
    if p.invert
        hcurr = plot([t0, tmax], (-med - currt) * [1 1], 'Color', [0.8 0.8 0.8]);
    else
        hcurr = plot([t0, tmax], (med + currt) * [1 1], 'Color', [0.8 0.8 0.8]);
    end
end
hold off
xlabel('Time (s)')
ylabel('Current (pA)')
title(sprintf('Spike detection: %s', fn));

% Prominence distribution
subplot(1,6,5)
hhist = histogram(proms, p.nbins);
[~, centers] = hist(proms, p.nbins);
xlabel('Prominence')
ylabel('Count')
title('Prominence threshold')

% Zoom in to spikes
subplot(1,6,6)
indsmat = ones(length(locs),1) * inds + locs * ones(1, length(inds));

% Remove 0s
if indsmat(1) < 0 
    oob = indsmat(:,1) < 1;
    indsmat(oob,:) = [];
end

% Remove extras
if indsmat(end) > l
    oob = indsmat(:,end) > l;
    indsmat(oob,:) = [];
end

wvs = d(indsmat);
hold on
if ~isempty(wvs)
    hwvs = plot(inds/fs, wvs', 'Color', [0.8 0.8 0.8]);
    hwvsmean = plot(inds/fs, mean(wvs,1), 'Color', 'k', 'LineWidth', 3);
end
plot([0 0], ylim(), 'k--')
hold off
xlabel('Time (ms)')
ylabel('Current (pA)')
title('Spikes')

%% UI elements
% Panel
hpan = uipanel(hfig,'Position',[0.0100 0.1000 0.0900 0.8500]);

% Display
toploc = [20 376 100 20];
uicontrol(hpan, 'Style', 'text', 'Position', toploc, 'String',...
    'Number of spikes:')
hsc = uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 16 0 0], 'String',...
    num2str(length(locs)));
uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 30 0 0], 'String',...
    'Firing rate:')
hfr = uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 46 0 0], 'String',...
    sprintf('%0.2f Hz', length(locs)/tmax));
uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 60 0 0], 'String',...
    'Half-max width:')
hws = uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 76 0 0], 'String',...
    sprintf('%0.2f ms', mean(ws)/fs));


% Prominence
uicontrol(hpan, 'Style', 'text', 'Position', [20 266 100 20], 'String', 'Prominence:')
uicontrol(hpan,'Style','edit','Position', [20 250 100 20], 'String',...
    num2str(promt), 'callback', @promchange);
function promchange(src, ~)
    % Prom
    promt = str2double(src.String);
    
    % New peaks
    [~, locs, ws, proms] = findpeaks(d2, 'MinPeakProminence', promt,...
    'MinPeakHeight', currt, 'MinPeakDistance', dis, 'MinPeakWidth', widtht);

    % Update counter
    hloc.XData = tvec(locs);
    hloc.YData = d(locs);
    
    % Update histograme
    hhist.BinCounts = hist(proms, centers);
    
    % Update means
    indsmat = ones(length(locs),1) * inds + locs * ones(1, length(inds));
    wvs = d(indsmat);
    subplot(1,6,6)
    delete(hwvs);
    delete(hwvsmean);
    hold on
    hwvs = plot(inds/fs, wvs', 'Color', [0.8 0.8 0.8]);
    hwvsmean = plot(inds/fs, mean(wvs,1), 'Color', 'k', 'LineWidth', 3);
    hold off
    
    % Update spike count and rate
    hsc.String = num2str(length(locs));
    hfr.String = sprintf('%0.2f Hz', length(locs)/tmax); 
    hws.String = sprintf('%0.2f ms', mean(ws)/fs);
end

% Current
uicontrol(hpan, 'Style', 'text', 'Position', [20 216 100 20], 'String', 'Current:')
uicontrol(hpan,'Style','edit','Position', [20 200 100 20], 'String',...
    num2str(currt), 'callback', @currchange);

function currchange(src, ~)
    % Prom
    currt = str2double(src.String);
    
    % New peaks
    [~, locs, ws, proms] = findpeaks(d2, 'MinPeakProminence', promt,...
    'MinPeakHeight', currt, 'MinPeakDistance', dis, 'MinPeakWidth', widtht);
    
    % Update the thresholds
    if p.subtractmovingmedian
        if p.invert
            hcurr.YData = -movmed - currt;
        else
            hcurr.YData = movmed + currt;
        end
    else
        if p.invert
            hcurr.YData = (-med - currt) * [1 1];
        else
            hcurr.YData = (med + currt) * [1 1];
        end
    end
    
    % Update counter
    hloc.XData = tvec(locs);
    hloc.YData = d(locs);
    
    % Update histograme
    hhist.BinCounts = hist(proms, centers);
    
    % Update means
    indsmat = ones(length(locs),1) * inds + locs * ones(1, length(inds));
    wvs = d(indsmat);
    subplot(1,6,6)
    delete(hwvs);
    delete(hwvsmean);
    hold on
    hwvs = plot(inds/fs, wvs', 'Color', [0.8 0.8 0.8]);
    hwvsmean = plot(inds/fs, mean(wvs,1), 'Color', 'k', 'LineWidth', 3);
    hold off
    
    % Update spike count and rate
    hsc.String = num2str(length(locs));
    hfr.String = sprintf('%0.2f Hz', length(locs)/tmax);
    hws.String = sprintf('%0.2f ms', mean(ws)/fs);
end

% Width
uicontrol(hpan, 'Style', 'text', 'Position', [20 166 100 20], 'String', 'Width:')
uicontrol(hpan, 'Style','edit','Position', [20 150 100 20], 'String',...
    num2str(widtht), 'callback', @widthchange);

function widthchange(src, ~)
    % Width
    widtht = str2double(src.String);
    
    % New peaks
    [~, locs, ws, proms] = findpeaks(d2, 'MinPeakProminence', promt,...
    'MinPeakHeight', currt, 'MinPeakDistance', dis, 'MinPeakWidth', widtht);
    
    % Update counter
    hloc.XData = tvec(locs);
    hloc.YData = d(locs);
    
    % Update histograme
    hhist.BinCounts = hist(proms, centers);
    
    % Update means
    indsmat = ones(length(locs),1) * inds + locs * ones(1, length(inds));
    wvs = d(indsmat);
    subplot(1,6,6)
    delete(hwvs);
    delete(hwvsmean);
    hold on
    hwvs = plot(inds/fs, wvs', 'Color', [0.8 0.8 0.8]);
    hwvsmean = plot(inds/fs, mean(wvs,1), 'Color', 'k', 'LineWidth', 3);
    hold off
    
    % Update spike count and rate
    hsc.String = num2str(length(locs));
    hfr.String = sprintf('%0.2f Hz', length(locs)/tmax);
    hws.String = sprintf('%0.2f ms', mean(ws)/fs);
end

% Save
uicontrol(hpan,'Style','pushbutton','Position', [20 100 100 20], 'String',...
    'Save', 'callback', @savedata);
function savedata(src, ~)
    % Parameter files
    if newparams
    params = struct('locs', locs, 't0', t0, 'tmax', tmax, 'fs', fs, 'promt',...
        promt, 'currt', currt, 'widtht', widtht, 'datawin', datawin, 'proms',...
        proms, 'dataparams', dataparams, 'finderparams', p, 'ws', ws);
    else
        params.locs = locs;
        params.proms = proms;
        params.promt = promt;
        params.currt = currt;
        params.widtht = widtht;
        params.finderparams = p;
        params.ws = ws;
    end
    save(fn_out, '-struct', 'params', '-v7.3');
    src.String = 'Data saved';
end

% Quit
uicontrol(hpan,'Style','pushbutton','Position', [20 80 100 20], 'String',...
    'Quit', 'callback', @quit);
function quit(~, ~)
    close(hfig);
end
end