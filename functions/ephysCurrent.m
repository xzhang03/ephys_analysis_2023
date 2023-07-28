function outputmat = ephysCurrent(fpath, varargin)
%ephysCurrent analyzes the currents from EPSC/IPSC events
%   ephysCurrent(fpath, varargin)
%% Parse
if nargin < 2
    varargin = {};
    if nargin < 1
        fpath = '';
    end
end

% Debug
% fpath = 'E:\ephys\stephen\May 2021 loose AVPV TH Stephen\211101\211101a AVPV THcre\211101a AVPV THcre_0002_preprocess.mat';

p = inputParser;

% Data handling parameters
addOptional(p, 'defaultpath', '\\anastasia\data\ephys\stephen\*_preprocess.mat');
addOptional(p, 'useprevparams', true); % Use previous parameters

% Analysis parameters
addOptional(p, 'channel', 1);
addOptional(p, 'EPSC', true);
addOptional(p, 'window', [-200 200]); % In ms
addOptional(p, 'bwin', [-100, -5]); % In ms, baseline window
addOptional(p, 'pwin', [0 100]); % In ms, peak window
addOptional(p, 'movingmedianwin', 20000) % Number of points in the moving median window

% Display parameters
addOptional(p, 'smoothwin', 5);
addOptional(p, 'pos', [50 200 1800 500]); % Figure position

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% IO
% Load
if isempty(fpath)
    [fn, fp] = uigetfile(p.defaultpath);
    [~, fn, ext] = fileparts(fn);
    
    % Full filename
    fpath = fullfile(fp, sprintf('%s%s', fn, ext));
else
    [fp, fn, ~] = fileparts(fpath);
end

% abf filename
fn_abf = fn(1:strfind(fn, '_preprocess')-1);
fpath_abf = fullfile(fp, sprintf('%s%s', fn_abf, '.abf'));

% Load parameters
params = load(fpath);
t0 = params.t0;
tmax = params.tmax;
fs = params.fs;
locs = params.locs;

% Load abf
[d, ~, ~] = abfload(fpath_abf);
d = d(params.datawin(1):params.datawin(2), p.channel);
movmed = movmedian(d, p.movingmedianwin);
l = length(d);

%% Indices
% Inds for EPSCs
inds = p.window(1)*fs : p.window(2)*fs;
inds_base = p.bwin(1)*fs : p.bwin(2)*fs;
inds_peak = p.pwin(1)*fs : p.pwin(2)*fs;
indsmat = ones(length(locs),1) * inds + locs * ones(1, length(inds));
indsmat_base = ones(length(locs),1) * inds_base + locs * ones(1, length(inds_base));
indsmat_peak = ones(length(locs),1) * inds_peak + locs * ones(1, length(inds_peak));

% Remove 0s
if indsmat(1) < 0
    oob = indsmat(:,1) <= 0;
    indsmat = indsmat(~oob, :);
    indsmat_base = indsmat_base(~oob, :);
    indsmat_peak = indsmat_peak(~oob, :);
    locs = locs(~oob);
end

% Remove extras
if indsmat(end,end) > l
    oob = indsmat(:,end) > l;
    indsmat = indsmat(~oob, :);
    indsmat_base = indsmat_base(~oob, :);
    indsmat_peak = indsmat_peak(~oob, :);
    locs = locs(~oob);
end

wvs = d(indsmat);
wvs_base = movmed(indsmat_base);
wvs_peak = d(indsmat_peak);

%% Current
% Current
tvec = t0 : 1/fs/1000 : tmax;
if p.EPSC
    amps = min(wvs_peak, [], 2) - mean(wvs_base,2);
else
    amps = max(wvs_peak, [], 2) - mean(wvs_base,2);
end
amps_smooth = smooth(amps, p.smoothwin);
amps_total =  mean(amps);

%% Stim
if isfield(params, 'stimmat')
    stimmat = params.stimmat;
    nstim = size(stimmat,1);
    dostim = true;
else
    dostim = false;
end

%% Plot
% Figure
figure('Position', p.pos);
% subplot(1,6,1:5)

hold on
plot(tvec(locs), amps_smooth)
plot([t0 tmax], amps_total * [1 1], 'k--')
if dostim
    for i = 1 : nstim
        plot([stimmat(i,1); stimmat(i,1)]/fs/1000, ylim()', '-r');
        plot([stimmat(i,1)+stimmat(i,2); stimmat(i,1)+stimmat(i,2)]/fs/1000, ylim()', '-r');
    end
end

hold off
xlabel('Time (s)')
ylabel('Intantaneous amplitude (pA)')
title(sprintf('Mean amplitude: %0.2f pA', amps_total))

% Zoom
%{
subplot(1,6,6)
hold on
if ~isempty(wvs)
    plot(inds/fs, wvs', 'Color', [0.8 0.8 0.8]);
    plot(inds/fs, mean(wvs,1), 'Color', 'k', 'LineWidth', 3);
end
plot([0 0], ylim(), 'k--')
hold off
xlabel('Time (ms)')
ylabel('Current (pA)')
title('EPSCs')
%}

%% Output
outputmat = zeros(601, size(stimmat,1));
for i = 1 : size(stimmat,1)
    t0 = stimmat(i,1) / 1000 / fs;
    trigt = (t0 - 60) : 0.2 : (t0 + 60);
    vq = interp1(tvec(locs), amps_smooth, trigt);
    outputmat(:,i) = vq;
end

end