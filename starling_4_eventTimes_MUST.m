% in this code I want to calculate the spaceRT, arrowRT and ITI and compare
% them with the ones I have from js psych files.

% trialStart: 1 , by fixation point: 500ms
% cardShow: 2 , 1000ms
% instructionMessage: 3 , no delay
% flipSpace: 4 , no delay <3000ms
% choiceAndFeedback: 5 , feedback is shown for 2000ms
% totalReward: 6 , shown for 1000ms

clc;
clear;
close all;

%%
totalRewardShowTime = 1000;
cardsShowTime = 1000;
choiceTime = 3500;
timeoutMessageTime = 2000;

%% Main folder
input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');
d = dir(input_folder);

isub = [d(:).isdir]; 
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = string(subFolders);

%%
% ptNumber = 1;

for p = 1:numel(ptIDs)
% for p = 8:8
    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);
    input_folder_pt = fullfile(input_folder, ptID);  

    % reading task data and extract ITI:
    bhvFiles = dir(fullfile(input_folder_pt, 'task_data*.csv'));
    bhvFile = fullfile(bhvFiles(1).folder, bhvFiles(1).name);
    bhvData = readtable(bhvFile);
    bhv_ITI = bhvData.interTrialInterval;
    timeout_idx = find(strcmp(bhvData.trialType, 'timeout'));

    % reading trigs data and calculate ITI:
    trigFile = fullfile(input_folder_pt, 'trigs.mat');
    load(trigFile, 'trigs');

    % Extract event timings
    trialStartTime         = trigs(2, trigs(1,:) == 1);
    cardShowTime           = trigs(2, trigs(1,:) == 2);
    instructionMessageTime = trigs(2, trigs(1,:) == 3);
    flipSpaceTime          = trigs(2, trigs(1,:) == 4);
    choiceAndFeedbackTime  = trigs(2, trigs(1,:) == 5);
    totalRewardTime        = trigs(2, trigs(1,:) == 6);
    
    % Make copies so we don’t overwrite directly
    cfTime = choiceAndFeedbackTime;
    trTime = totalRewardTime;
    
    % Insert NaN and shift for timeout indices
    for k = 1:numel(timeout_idx)
        idx = timeout_idx(k);
    
        % Choice & Feedback
        cfTime = [cfTime(1:idx-1), NaN, cfTime(idx:end)];
    
        % Total Reward
        trTime = [trTime(1:idx-1), NaN, trTime(idx:end)];
    end

    choiceAndFeedbackTime = cfTime;
    totalRewardTime = trTime;

    nTrials = length(trialStartTime);
    photodiodeITI = nan(1, nTrials);

    for i = 1:nTrials-1
        if ismember(i, timeout_idx)
            trialEndTime = flipSpaceTime(i) + choiceTime + timeoutMessageTime + cardsShowTime;
        else
            trialEndTime = totalRewardTime(i) + totalRewardShowTime ;
        end
        photodiodeITI(i+1) = trialStartTime(i+1) - trialEndTime;
    end

    save(fullfile(input_folder_pt, 'eventTimes.mat'), ...
        'trialStartTime', ...
        'cardShowTime', ...
        'instructionMessageTime', ...
        'flipSpaceTime', ...
        'choiceAndFeedbackTime', ...
        'totalRewardTime', ...
        'photodiodeITI');

end
