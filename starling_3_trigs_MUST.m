% in this code I'm trying to make a trigs mat files for each 
% patient regarding both "response" and "timeout" trials.
% there are each trig code:
% trialStart: 1 , by fixation point: 500ms
% cardShow: 2 , 1000ms
% instructionMessage: 3 , no delay
% flipSpace: 4 , no delay <3000ms
% choiceAndFeedback: 5 , feedback is shown for 2000ms
% totalReward: 6 , shown for 1000ms

clc;
clear;
close all;
% add new patients here
differentPatients = {'202421', '202509', '202511', '202512'};
%% Main folder
input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');
d = dir(input_folder);


isub = [d(:).isdir]; 
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = string(subFolders);

%%
% for p = 1:numel(ptIDs)
for p = 7:7
    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);

    output_folder = fullfile(input_folder, ptID);

    % --- find NEV file ---
    if any(strcmp(ptID, differentPatients ))
        nevList = dir(fullfile(input_folder_pt, '*.nev'));
    else
        nevList = dir(fullfile(input_folder_pt, '\hub_neural_data\*.nev'));
    end
    
    if length(nevList) > 1
        error('many nev files available for this patient. Please specify...')
    elseif isempty(nevList)
        error('no nev files found...')
    else
        nevFile = fullfile(nevList.folder, nevList.name);
    end
    
    %  get photodiode 
    [nevPath,nevName,~] = fileparts(nevFile);
    NS5 = openNSx(fullfile(nevPath,[nevName '.ns5']));
    original_freq = NS5.MetaTags.SamplingFreq;

    desired_freq = 1000; % Hz
    downsample_steps = original_freq/desired_freq;
    lp_cutoff = 60; % low-pass filter 60 Hz
    [b,a] = butter(4, lp_cutoff/(original_freq/2), 'low');
    photodiode = double(NS5.Data(2, :));
    photodiode = filtfilt(b,a, photodiode);
    % downsampling
    % photodiode = photodiode(1:downsample_steps:end); Elliot recommended
    % the resample function
    
    [p, q] = rat(desired_freq / original_freq);
    photodiode = resample(photodiode, p, q);

    
    %  rising edges
    thresh_diff = 500; % relative threshold
    window_ms = 20;    % lookback window
    
    rising_idx = [];
    for i = (window_ms+1):length(photodiode)
        baseline = mean(photodiode(i-window_ms:i-1));
        if (photodiode(i) - baseline) > thresh_diff
            if isempty(rising_idx) || i - rising_idx(end) > window_ms
                rising_idx(end+1) = i;
            end
        end
    end
    
    % falling edges corresponding to each rising edge 

    
    falling_idx = nan(size(rising_idx));
    for r = 1:length(rising_idx)
        start_i = rising_idx(r);
        % search forward until signal drops back near baseline
        base = mean(photodiode(max(1,start_i-window_ms):start_i));
        drop_thresh = base + thresh_diff/4;
        j = start_i+1;
        while j <= length(photodiode) && photodiode(j) > drop_thresh
            j = j+1;
        end
        if j <= length(photodiode)
            falling_idx(r) = j;
        end
    end
    
    % remove short pulses <50 ms 
    fallRiseBadGap = 50;
    valid_mask = (falling_idx - rising_idx) >= fallRiseBadGap;
    rising_idx = rising_idx(valid_mask);
    falling_idx = falling_idx(valid_mask);
    
    % detect doublets 
    isi = diff(rising_idx);
    doublet_idx = find(50 < isi & isi <= 200);
    doublet_rising_idx = rising_idx(doublet_idx);
    
    % count rising edges between doublets
    n_between_list = nan(1, length(doublet_idx)-1);
    for k = 1:length(doublet_idx)-1
        end_curr   = rising_idx(doublet_idx(k)+1);
        start_next = rising_idx(doublet_idx(k+1));
        between_edges = rising_idx(rising_idx > end_curr & rising_idx < start_next);
        n_between_list(k) = numel(between_edges);
    end
    
    % remove bad doublets
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
    
    % spikes between doublets 
    n_between_doublets = nan(1, length(doublet_rising_idx));
    spike_times_between = cell(1, length(doublet_rising_idx));
    
    for k = 1:length(doublet_rising_idx)
        pos = find(rising_idx == doublet_rising_idx(k), 1, 'first');
        end_curr = rising_idx(pos+1);
        
        if k < length(doublet_rising_idx)
            start_next = doublet_rising_idx(k+1);
        else
            start_next = length(photodiode);
        end
        
        between_edges = rising_idx(rising_idx > end_curr & rising_idx < start_next);
        n_between_doublets(k) = numel(between_edges);
        spike_times_between{k} = between_edges;
    end
    
    %  trigs structure 
    all_trigs = [];  % 2×all trigs first is trig, second is time in ms
    
    for k = 1:length(doublet_rising_idx)
        this_trigs = [];
        this_times = [];
        
        this_trigs(end+1) = 1;
        this_times(end+1) = double(doublet_rising_idx(k)); 
        
        if n_between_doublets(k) >= 5
            % full response trial
            codes = [2 3 4 5 6];
            times = double(spike_times_between{k}(1:5));
            
        elseif n_between_doublets(k) == 4
            % overlapped trial
            codes = [2 3 4 5 6];
            t = double(spike_times_between{k});
            times = [t(1:2), t(2)+30, t(3:4)];
            
        elseif n_between_doublets(k) == 3
            % 3 spike trial
            codes = [2 3 4];
            t = double(spike_times_between{k});
            times = t(1:3);

            
        elseif n_between_doublets(k) == 2
            % timeout trial
            codes = [2 3 4];
            t = double(spike_times_between{k});
            times = [t(1), t(2), t(2)+30];
        else
            codes = [];
            times = [];
        end
        
        this_trigs = [this_trigs codes];
        this_times = [this_times times];
        
        all_trigs = [all_trigs [this_trigs; this_times]];
    end
    
    trigs = all_trigs;  
    save(fullfile(output_folder, 'trigs.mat'), 'trigs');
    
end



%%

    % --- handle photodiode error for ptID 202510 ---
    % if strcmp(ptID, '202510')
    %     photodiode(681400:690553) = 0;
    %     photodiode(996443:1021620) = 0;
    %     photodiode(682984:1021370) = -photodiode(682984:1021370);
    %     photodiode(682984:1021370) = photodiode(682984:1021370) * 25;
    % end
    % 