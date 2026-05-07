% in this code I'm trying to find out whether the trials I've found align
% with the task_data file
% result: 202510 trial 63 does not match it's because of the noise on data
% n_between_doublets are the number of spikes of photodiodes between each
% doublet found
clc;
clear;
close all;

% add new patients here
differentPatients = {'202514', '202518', '202601'};
%%
input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');
d = dir(input_folder);

isub = [d(:).isdir]; 
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = subFolders; 
ptIDs = string(ptIDs);

%% 
for p = 1:numel(ptIDs)
% for p = 9:9
    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);

    output_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw', ptID);
    input_bhv_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw', ptID);
    
    
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
    
    % load data from ns5 to get photodiode:
    [nevPath,nevName,nevExt] = fileparts(nevFile);
    NS5 = openNSx(fullfile(nevPath,[nevName '.ns5']));
    original_freq = NS5.MetaTags.SamplingFreq;
    desired_freq = 1000;

    % photodiode
    downsample_steps = original_freq/desired_freq;
    photodiode = double(NS5.Data(2, 1:downsample_steps:end));
    lp_cutoff = 60;
    [b,a] = butter(4, lp_cutoff/(desired_freq/2), 'low');
    photodiode = filtfilt(b,a, photodiode);
    
    % handling photodiode error in 202510
    if strcmp(ptID, '202510')
        photodiode(681400:695553) = 0;
        photodiode(996443:1021620) = 0;
        photodiode(682984:1021370) = -photodiode(682984:1021370);
        photodiode(682984:1021370) = photodiode(682984:1021370) * 25;
    end
    
    
    
    % detect rising edges & doublets
    thresh_diff = 500;              % threshold relative to baseline
    window_ms = 20;                 % lookback window
    
    rising_idx = [];                % store detected rising edges
    
    for i = (window_ms+1):length(photodiode)
        baseline = mean(photodiode(i-window_ms:i-1));  % average in 20 ms before
        if (photodiode(i) - baseline) > thresh_diff
            % rising edge if difference exceeds +500
            if isempty(rising_idx) || i - rising_idx(end) > window_ms
                % prevent multiple detections within the same edge
                rising_idx(end+1) = i;
            end
        end
    end
    
    
    isi = diff(rising_idx);
    doublet_idx = find(50 < isi & isi <= 200); % finding doublets as the begining of trials
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
    
    
    % figure;
    % plot(photodiode, 'k');  
    % hold on;
    % 
    % plot(doublet_rising_idx, photodiode(doublet_rising_idx), 'ro', ...
    %      'MarkerSize', 6, 'LineWidth', 1.5);
    
    % count rising edges between each doublet_rising_idx for finding normal/timeout trials
    n_between_doublets = nan(1, length(doublet_rising_idx));
    
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
    end
    
    
    % build trial type labels
    trialType = repmat("response", length(n_between_doublets), 1); % default
    trialType(n_between_doublets == 2 | n_between_doublets == 3) = "timeout";
    
    T = table((1:length(n_between_doublets))', n_between_doublets(:), trialType, ...
              'VariableNames', {'Trial', 'n_between_doublets', 'trialType'});
    
    csvFile = fullfile(output_folder, 'trialsInfo.csv');
    writetable(T, csvFile);
    % comparing trials for all the data
    % find csv file starting with 'task_data'
    bhvFiles = dir(fullfile(input_bhv_folder, 'task_data*.csv'));
    if isempty(bhvFiles)
        error('No task_data file found in %s', input_bhv_folder);
    end
    bhvFile = fullfile(bhvFiles(1).folder, bhvFiles(1).name);
    
    % read behavioral file
    bhvData = readtable(bhvFile);
    
    % copy trialType column into new variable
    bhv_trialType_raw = bhvData.trialType;
   
    
    % add raw trial types into your results table
    T.trialType_raw = strings(height(T),1);
    nToCompare = min(length(bhv_trialType_raw), height(T));
    T.trialType_raw(1:nToCompare) = string(bhv_trialType_raw(1:nToCompare));
    
    % write updated table back
    csvFile = fullfile(output_folder, 'trialsInfo.csv');
    writetable(T, csvFile);
    
    % compare trial types
    if all(T.trialType(1:nToCompare) == T.trialType_raw(1:nToCompare))
        fprintf('Trial types match\n');
    else
        fprintf('Trial types do not match\n');
    end
end