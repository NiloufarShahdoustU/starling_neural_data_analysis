clc;
clear;
close all;
ptID = '202601'; 
% add pts with new daq system
differentPatients = {'202514', '202518', '202601'};
%% 
if any(strcmp(ptID, differentPatients ))
    nevList = dir(fullfile('\\155.100.91.44\d\Data\Nill\starling\raw', ptID, '\nsp_photodiode_data\*.nev'));
else
    nevList = dir(fullfile('\\155.100.91.44\d\Data\Nill\starling\raw', ptID, '*.nev'));
end

if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end

%% load data from ns5 to get photodiode:
[nevPath,nevName,nevExt] = fileparts(nevFile);
NS5 = openNSx(fullfile(nevPath,[nevName '.ns5']));
original_freq = NS5.MetaTags.SamplingFreq;
%% load neural data from ns2:
% NS2 = openNSx(fullfile(nevPath,[nevName '.ns2']));
% desired_freq = NS2.MetaTags.SamplingFreq;

desired_freq = 1000; % Hz
downsample_steps = original_freq/desired_freq;
lp_cutoff = 30; % low-pass filter 30 Hz
[b,a] = butter(4, lp_cutoff/(desired_freq/2), 'low');
photodiode = double(NS5.Data(2, :));
photodiode = filtfilt(b,a, photodiode);
% downsampling
photodiode = photodiode(1:downsample_steps:end);

%% handling photodiode error in 202510
if strcmp(ptID, '202510')
    photodiode(681400:695553) = 0;
    photodiode(996443:1021620) = 0;
    photodiode(682984:1021370) = -photodiode(682984:1021370);
    photodiode(682984:1021370) = photodiode(682984:1021370) * 25;
end



%% detect rising edges & doublets
thresh_diff = 500;              % threshold relative to baseline
window_ms = 20;                 % lookback window (ms)
fs = desired_freq;              % sampling rate (Hz)
window_samples = round(window_ms * 1e-3 * fs);

rising_idx = [];                % store detected rising edges

for i = (window_samples+1):length(photodiode)
    baseline = mean(photodiode(i-window_samples:i-1));  % average in 30 ms before
    if (photodiode(i) - baseline) > thresh_diff
        % rising edge if difference exceeds +500
        if isempty(rising_idx) || i - rising_idx(end) > window_samples
            % prevent multiple detections within the same edge
            rising_idx(end+1) = i;
        end
    end
end


isi = diff(rising_idx);
doublet_idx = find(50 < isi & isi <= 200);
doublet_rising_idx = rising_idx(doublet_idx);

% count rising edges between doublets
n_between_list = nan(1, length(doublet_idx)-1);

for k = 1:length(doublet_idx)-1
    % define edges to exclude:
    end_curr   = rising_idx(doublet_idx(k)+1);   % 2nd rising of current doublet
    start_next = rising_idx(doublet_idx(k+1));   % 1st rising of next doublet
    
    % rising edges strictly between them
    between_edges = rising_idx(rising_idx > end_curr & rising_idx < start_next);
    n_between_list(k) = numel(between_edges);
end


    % --- detect falling edges corresponding to each rising edge ---
    falling_idx = nan(size(rising_idx));
    for r = 1:length(rising_idx)
        start_i = rising_idx(r);
        % search forward until signal drops back near baseline
        base = mean(photodiode(max(1,start_i-50):start_i));
        drop_thresh = base + thresh_diff/4; % more lenient threshold for falling
        j = start_i+1;
        while j <= length(photodiode) && photodiode(j) > drop_thresh
            j = j+1;
        end
        if j <= length(photodiode)
            falling_idx(r) = j;
        end
    end
    
    % --- remove short pulses (<50 ms = 50 samples at 1000 Hz) ---
    valid_mask = (falling_idx - rising_idx) >= 50;
    rising_idx = rising_idx(valid_mask);
    falling_idx = falling_idx(valid_mask);

%% 
% Find indices of bad doublets
trials_good = 5;
trials_good_miss = 3;
trials_message_show_and_space_overlaped = 4;
trials_message_show_and_space_overlaped_miss = 2;
consecutive = 0;

bad_idx = find(n_between_list ~= trials_good & ...
               n_between_list ~= trials_good_miss & ...
               n_between_list ~= trials_message_show_and_space_overlaped & ...
               n_between_list ~= trials_message_show_and_space_overlaped_miss & ...
               n_between_list ~= consecutive);

doublet_rising_idx(bad_idx+1) = [];


%%
figure;
plot(photodiode, 'k');  
hold on;

plot(doublet_rising_idx, photodiode(doublet_rising_idx), 'ro', ...
     'MarkerSize', 6, 'LineWidth', 1.5);

% Add trial numbers on top of each red circle
for t = 1:length(doublet_rising_idx)
    text(doublet_rising_idx(t), photodiode(doublet_rising_idx(t)) + 200, ... % +200 shifts text above the circle
         num2str(t), ...
         'HorizontalAlignment', 'center', ...
         'FontSize', 8, ...
         'Color', 'b', ...
         'FontWeight', 'bold');
end

