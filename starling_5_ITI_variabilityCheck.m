clc;
clear;
close all;
%% Main folder
input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');
d = dir(input_folder);

isub = [d(:).isdir]; 
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = string(subFolders);

%% Prepare figure
nPatients = numel(ptIDs);
nCols = ceil(sqrt(nPatients));
nRows = ceil(nPatients / nCols);

figure('Position', [100 100 1600 900]);

for p = 1:nPatients
    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);
    input_folder_pt = fullfile(input_folder, ptID);  

    % reading task data and extract ITI:
    bhvFiles = dir(fullfile(input_folder_pt, 'task_data*.csv'));
    if isempty(bhvFiles)
        warning('No behavioral file found for %s. Skipping.', ptID);
        continue;
    end
    bhvFile = fullfile(bhvFiles(1).folder, bhvFiles(1).name);
    bhvData = readtable(bhvFile);
    bhv_ITI = bhvData.interTrialInterval;

    % reading eventTimes:
    eventTimesFile = fullfile(input_folder_pt, 'eventTimes.mat');
    eventTimes = load(eventTimesFile);

    ITI_variability = eventTimes.photodiodeITI - bhv_ITI';

    medVal = median(ITI_variability, 'omitnan');
    madVal = mad(ITI_variability, 1); 
    threshold = 3; 
    keepIdx = abs(ITI_variability - medVal) <= threshold * madVal;
    ITI_clean = ITI_variability(keepIdx);

    fprintf('min of variability: %.2f\n', min(ITI_clean));
    fprintf('max of variability: %.2f\n', max(ITI_clean));
    fprintf('mean of variability: %.2f\n', mean(abs(ITI_clean)));
    fprintf('std of variability: %.2f\n', std(ITI_clean));
    fprintf('mode of variability: %.4f\n', mode(abs(ITI_clean)));
    fprintf('median of variability: %.4f\n', median(abs(ITI_clean)));

    subplot(nRows, nCols, p);
    plot(ITI_clean);
    title(ptID, 'Interpreter', 'none');

end

