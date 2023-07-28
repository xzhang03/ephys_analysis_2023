function outputmat = ephysInstaFR(fpath, varargin)
% ephysInstaFR plots the instattaneous firing rate
% ephysInstaFR(fpath, varargin)

%% Parse
if nargin < 2
    varargin = {};
    if nargin < 1
        fpath = '';
    end
end

% Debug
% fpath = 'E:\ephys\stephen\May 2021 loose AVPV TH Stephen\210513\210513a loose AVPV TH Stephen\210513a loose AVPV TH Stephen_0001_preprocess.mat';

p = inputParser;

% Data handling parameters
addOptional(p, 'defaultpath', '\\anastasia\data\ephys\stephen\*_preprocess.mat');
addOptional(p, 'useprevparams', true); % Use previous parameters

% Analysis parameters


% Display parameters
addOptional(p, 'smoothwin', 5);
addOptional(p, 'pos', [50 200 1800 500]); % Figure position
addOptional(p, 'ylim', [0 6]);

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
    
    % Full filename
    fpath = fullfile(fp, fn);
end

% Load
params = load(fpath, '-mat');

% Read parameters
if p.useprevparams && isfield(params, 'smoothwin')
    smoothwin = params.smoothwin;
    disp('Using previous parameters');
    newparams = false;
else
    smoothwin = p.smoothwin;
    newparams = true;
end

%% Spikes
% Spike info
t0 = params.t0;
tmax = params.tmax;
fs = params.fs;
locs = params.locs;

% FR
tvec = t0 : 1/fs/1000 : tmax;
fr = 1./smooth(diff(tvec(locs)), smoothwin);
fr_total = length(locs)/tmax;

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

hold on
plot(tvec(locs(2:end)), fr)
plot([t0 tmax], fr_total * [1 1], 'k--')
if dostim
    for i = 1 : nstim
        plot([stimmat(i,1); stimmat(i,1)]/fs/1000, p.ylim', '-r');
        plot([stimmat(i,1)+stimmat(i,2); stimmat(i,1)+stimmat(i,2)]/fs/1000, p.ylim', '-r');
    end
end

hold off
ylim(p.ylim)
xlabel('Time (s)')
ylabel('Intantaneous firing rate (Hz)')
title(sprintf('Firing rate: %0.2f Hz', fr_total))

%% Save
if newparams
    params.smoothwin = smoothwin;
    save(fpath, '-struct', 'params', '-v7.3');
end

%% Output
outputmat = zeros(601, size(stimmat,1));
for i = 1 : size(stimmat,1)
    t0 = stimmat(i,1) / 1000 / fs;
    trigt = (t0 - 60) : 0.2 : (t0 + 60);
    vq = interp1(tvec(locs(2:end)), fr, trigt);
    outputmat(:,i) = vq;
end

end