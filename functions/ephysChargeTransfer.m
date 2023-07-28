function ephysChargeTransfer(fpath, varargin)
%ephysChargeTransfer analyzes the charge transfers from EPSC/IPSC events
%   ephysChargeTransfer(fpath, varargin)

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
addOptional(p, 'window', [-30 300]); % In ms
addOptional(p, 'risewindow', 8); % Window for the rise phase of IPSC/EPSCs in ms
addOptional(p, 'CTpercentile', 0.01); % 0.01 means the window is the width of the 1% max values
addOptional(p, 'smoothwin', 50); % Smooth window, in points. Set to 0 if no smoothing.

% Display parameters
addOptional(p, 'pos', [50 200 1800 500]); % Figure position
addOptional(p, 'nbins', 20);
addOptional(p, 'checkresults', false); % Figure to check results

% Save
addOptional(p, 'save', true);

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

% Loca
locs = params.locs + params.datawin(1) - 1;
n = length(locs);

% Load abf
[d, si, ~] = abfload(fpath_abf);
d = d(:, p.channel);
ppms = 1000/si; % pts/ms

% Smooth
if p.smoothwin > 0
    d = smooth(d, p.smoothwin);
end

%% Basic numbers
% Point window
windowpt = p.window * ppms;
risewindowpt = p.risewindow * ppms;

%% Datamatrix
% Indicies
inds1 = ones(n,1) * (windowpt(1) : windowpt(2));
inds2 = locs * ones(1, windowpt(2) - windowpt(1) + 1);
inds = inds1 + inds2;

% Remove out of bound on the left side
if inds(1,1) <= 0
    oob = inds(:,1) <= 0;
    inds = inds(~oob, :);
    locs = locs(~oob);
    n = size(inds,1);
end

% Remove out of bound
if inds(end,end) > length(d)
    oob = inds(:,end) > length(d);
    inds = inds(~oob, :);
    locs = locs(~oob);
    n = size(inds,1);
end


% Matrix
d2 = d(inds');

% Remove points that belong to the next EPSC/IPSC
for i = 1 : n
    locstmp = locs - risewindowpt;
    locstmp(i) = 0;
    if any(ismember(locstmp, inds(i,:)))
        indsintercept = find(ismember(locstmp, inds(i,:)));
        for j = 1 : length(indsintercept)
            indsintercept_curr = indsintercept(j);
            badind = find((locs(indsintercept_curr)-risewindowpt) == inds(i,:));
            d2(badind:end, i) = nan;
        end
    end
end

% Resting current
yrst = nanmean(d2(1: -windowpt(1)-risewindowpt,:));

%% Calculations
if p.useprevparams && isfield(params, 'CT') && isfield(params, 'CTwin') && isfield(params, 'Imax')
    disp('Using previous parameters.')
    CT = params.CT;
    CTwin = params.CTwin;
    Imax = params.Imax;
    p.save = false;
else
    % Imax
    Imax = (d2(-windowpt(1),:) - yrst)';

    % Charge transfers
    CT = zeros(n, 1);
    CTwin = zeros(n, 2);

    for i = 1 : n
        % Threshold
        thresh = Imax(i) * p.CTpercentile + yrst(i);

        % Left index
        ind1 = find(d2(1:-windowpt(1),i) < thresh, 1, 'last');
        if isempty(ind1)
            [~, ind1] = min(d2(1:-windowpt(1),i));
        end

        % Right index
        ind2 = find(d2(-windowpt(1):end,i) < thresh, 1, 'first');
        if isempty(ind2)
            [~, ind2] = min(d2(-windowpt(1):end,i));
        end
        ind2 = ind2 - windowpt(1) - 1;

        % Calculate
        CTwin(i,1) = ind1;
        CTwin(i,2) = ind2;
        CT(i) = sum(d2(ind1:ind2, i) - yrst(i));


    %     % Debug
    %     plot(d2(:,1))
    %     hold on
    %     plot([ind1 ind2], d2([ind1 ind2],1), 'r.')
    %     hold off
    end

    % Charge transfer
    CT = CT / ppms;
end
%% Plot
figure('Position', p.pos);

% Traces
subplot(1,6,1:4)
plot(d2, 'Color', [0.8 0.8 0.8])
hold on
plot(nanmean(d2,2), 'Color', [0 0 0], 'LineWidth', 2);
plot(xlim(), nanmean(yrst) * [1 1], 'k-');
hold off
xlabel('Time')
ylabel('Current (pA)')

% Currents
subplot(1,6,5)
histogram(Imax, p.nbins);
xlabel('Imax (pA)')
ylabel('Count')
title('Current')

% Charge transfers
subplot(1,6,6)
histogram(CT, p.nbins);
xlabel('Charge transfer (pA * ms)')
ylabel('Count')
title('Charge transfer')

% display
fprintf('Mean Imax (pA): %0.2f.\n', nanmean(Imax))
fprintf('Mean Charge transfered (pA * ms): %0.2f.\n', nanmean(CT))

%% Check result
if p.checkresults
    figure

    % Normalize
    d2n = d2;
    for i = 1 : n
        d2n(:,i) = mat2gray(d2(:,i)) + i;
        d2n(d2n(:,i) == i+1, i) = nan;
    end

    % 4 Quartiles
    nq = ceil(n/4);

    for i = 1 : 4
        % Subplot
        subplot(1,4,i);

        % Inds
        istart = (i-1) * nq + 1;
        if i == 4
            iend = n;
        else
            iend = i * nq;
        end

        plot(d2n(:,istart:iend) - istart, 'Color', [0.5 0.5 0.5]);

        % Lims
        ylim([0 nq])
        xlim([0 windowpt(2)-windowpt(1)])

        if i > 1
            set(gca, 'YTickLabel', {})
        end

        % Mark ponits
        hold on
        for j = istart : iend
            plot(CTwin(j,:), d2n(CTwin(j,:),j) - istart,'r-')
        end
        hold off
    end
end

%% Save
if p.save
    params.CT = CT;
    params.CTwin = CTwin;
    params.Imax = Imax;
    
    save(fpath, '-struct', 'params', '-v7.3');
    disp('Data saved');
end


end

