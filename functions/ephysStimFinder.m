function fn_out = ephysStimFinder(fpath, varargin)
% ephysLoosePatchSpikeFinder finds stim time points in loose patch data
% fn_out = ephysStimFinder(fpath, varargin)

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
addOptional(p, 'usedatawindow', true); % Use a window of data
addOptional(p, 'channel', 3);
addOptional(p, 'threshold', 1000);

% Display parameters
addOptional(p, 'pos', [50 200 1800 500]); % Figure position

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
    de = d(datawin(1) : datawin(2),1);
    d = d(datawin(1) : datawin(2),p.channel);
    
    close(hdatawin)
elseif p.usedatawindow && ~newparams
    datawin = params.datawin;
    de = d(datawin(1) : datawin(2),1);
    d = d(datawin(1) : datawin(2),p.channel);
else
    datawin = [1 length(d)];
    de = d(:,1);
    d = d(:,p.channel);
end

%% Preprocess
% Get channel
l = length(d);

% Sampling rate and times
fs = 1000/si; % in ms-1
tmax = l/fs/1000;
t0 = 1/fs/1000;
tvec = t0 : 1/fs/1000 : tmax;

% Onsets and offsets
stimmat = chainfinder(d > p.threshold);
if isempty(stimmat)
    disp('No stim found.');
    return;
end
ons = stimmat(:,1);

%% Figure
hfig = figure('Position', p.pos);

subplot(2,1,1)
hold on
plot(tvec, d(:,1));
plot(tvec(ons), d(ons), 'r.');
hold off
ylabel('Current (pA)')
title(sprintf('Stim detection: %s', fn));

subplot(2,1,2)
hold on
plot(tvec, de);
plot(tvec(ons), de(ons), 'ro');
hold off
xlabel('Time (s)')
ylabel('Current (pA)')

%% UI elements
% Panel
hpan = uipanel(hfig,'Position',[0.0100 0.1000 0.0900 0.8500]);

% Display
toploc = [20 376 100 20];
uicontrol(hpan, 'Style', 'text', 'Position', toploc, 'String',...
    'Number of stims:')
hsc = uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 16 0 0], 'String',...
    num2str(size(stimmat,1)));
uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 30 0 0], 'String',...
    'Stim rate:')
hfr = uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 46 0 0], 'String',...
    sprintf('%0.2f Hz', size(stimmat,1)/tmax));
uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 60 0 0], 'String',...
    'Width:')
hws = uicontrol(hpan, 'Style', 'text', 'Position', toploc - [0 76 0 0], 'String',...
    sprintf('%0.2f ms', mean(stimmat(:,2))/fs));

% Save
uicontrol(hpan,'Style','pushbutton','Position', [20 100 100 20], 'String',...
    'Save', 'callback', @savedata);
function savedata(src, ~)
    % Parameter files
    if newparams
        % Fixed inthe future but really shouldn't be here
    else
        params.stimmat = stimmat;
        params.stimtime = stimmat / fs / 1000;
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