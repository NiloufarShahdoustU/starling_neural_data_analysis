% outcome based spectogram

% trialStart: 1 , by fixation point: 500ms
% cardShow: 2 , 1000ms
% instructionMessage: 3 , no delay
% flipSpace: 4 , no delay <3000ms
% choiceAndFeedback: 5 , feedback is shown for 2000ms
% totalReward: 6 , shown for 1000ms

% outcome based spectogram

clc;
clear;
close all;

differentPatients = {'202421', '202509', '202511', '202512'};
needsConcatenationPts = {'202514','202521'};

%% all these times are in ms
choiceFeedbackWindow = 1500;

%% Main folder
output_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\spectrograms\outcome\');
input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');

d = dir(input_folder);

isub = [d(:).isdir]; 
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = string(subFolders);

%% Loop through participants

for p = 1:numel(ptIDs)

    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);
    
    input_folder_pt = fullfile(input_folder, ptID); 

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % reading neural data

    if any(strcmp(ptID, differentPatients))
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
    
    % ns2 

    [nevPath, nevName, ~] = fileparts(nevFile);
    NS2 = openNSx(fullfile(nevPath, [nevName '.ns2']));
    original_freq = NS2.MetaTags.SamplingFreq;

    % reading selected channels using ptTrodesSTARLING

    [trodeLabels, isECoG, ~, ~, anatomicalLocs] = ptTrodesSTARLING(ptID);

    selectedChans = find(isECoG);
    selectedChans = selectedChans(1:end-1);
    
    % remove channels whose anatomical location contains white matter or NaC
    SelectedAnatomicalLoc = anatomicalLocs(selectedChans);
    
    badLoc = contains(SelectedAnatomicalLoc, 'white matter', 'IgnoreCase', true) | ...
             contains(SelectedAnatomicalLoc, 'NaC', 'IgnoreCase', true);
    
    selectedChans = selectedChans(~badLoc);
    SelectedAnatomicalLoc = anatomicalLocs(selectedChans);

    nChans = length(selectedChans);

    % Load and filter LFP data

    LFPData = [];

    [b1, a1] = iirnotch(60/(original_freq/2),  (60/(original_freq/2))/25);
    [b2, a2] = iirnotch(120/(original_freq/2), (120/(original_freq/2))/25);
    
    for ch = 1:nChans

        if ismember(ptID, needsConcatenationPts)

            tmpData = [];

            for c = 1:numel(NS2.Data)
                tmpData = [tmpData, double(NS2.Data{c}(selectedChans(ch), :))];
            end
    
            tmp = filtfilt(b1, a1, tmpData);
            tmp = filtfilt(b2, a2, tmp);

            LFPData(ch, :) = tmp;

            clear tmpData tmp
    
        else

            tmp = filtfilt(b1, a1, double(NS2.Data(selectedChans(ch), :)));
            tmp = filtfilt(b2, a2, tmp);

            LFPData(ch, :) = tmp;

            clear tmp

        end
    end

    % Denoise data using Elliot's function

    LFPData = remove1stPC(LFPData);

    % optional common average rereference
    % LFPData = LFPData - mean(LFPData, 2);

    % reading task data

    bhvFiles = dir(fullfile(input_folder_pt, 'task_data*.csv'));

    if isempty(bhvFiles)
        error('No behavioral file found for ptID %s', ptID);
    end

    bhvFile = fullfile(bhvFiles(1).folder, bhvFiles(1).name);
    bhvData = readtable(bhvFile);
    nTrials = size(bhvData, 1);

    % reading eventTimes

    eventTimesFile = fullfile(input_folder_pt, 'eventTimes.mat');

    if ~exist(eventTimesFile, 'file')
        error('No eventTimes.mat found for ptID %s', ptID);
    end

    eventTimes = load(eventTimesFile);

    choiceOutcome = eventTimes.choiceAndFeedbackTime;
    trialStart = eventTimes.trialStartTime;
    
    % Spectrogram parameters

    fWin = [1 200];
    waitBar = 0;
    motherWaveletParam = 6;

    Sbl_all = cell(1, nChans);
    freq_all = cell(1, nChans);
    
    % Spectrogram calculation

    for ch = 1:nChans

        fprintf('\nDoing spectral calculations for chan %d of %d', ch, nChans)
    
        clear S
        S = [];
        baselineSum = [];
        baselineCount = 0;
    
        for tt = 1:nTrials
    
            if ~isnan(choiceOutcome(tt))
    
                % main LFP segment around choice/feedback

                whichData = (choiceOutcome(tt) - choiceFeedbackWindow):(choiceOutcome(tt) + choiceFeedbackWindow - 1);

                if min(whichData) < 1 || max(whichData) > size(LFPData, 2)
                    warning('Skipping trial %d for ptID %s, channel %d: main window out of bounds.', tt, ptID, ch);
                    continue
                end

                LFPseg = LFPData(ch, whichData);
    
                [W, period, scale] = basewaveERP(LFPseg, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
                pow = single(10 * log10(abs(W).^2));
    
                S(:, :, tt) = pow;
    
                if isempty(freq_all{ch})
                    freq_all{ch} = 1 ./ period;
                end
    
                % baseline segment around trial start

                baseData = (trialStart(tt) - 750):(trialStart(tt) + 500 - 1);

                if min(baseData) < 1 || max(baseData) > size(LFPData, 2)
                    warning('Skipping baseline for trial %d for ptID %s, channel %d: baseline window out of bounds.', tt, ptID, ch);
                    continue
                end

                LFPbase = LFPData(ch, baseData);
    
                [bW, bperiod, bscale] = basewaveERP(LFPbase, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
                bPow = single(10 * log10(abs(bW).^2));
    
                baselineWindow = bPow(:, 100:end-150);
    
                if isempty(baselineSum)
                    baselineSum = zeros(size(bPow, 1), 1, 'single');
                end
    
                baselineSum = baselineSum + mean(baselineWindow, 2);
                baselineCount = baselineCount + 1;

            end
        end
    
        % baseline correction

        if baselineCount == 0
            warning('No valid baseline trials for ptID %s, channel %d. Saving empty Sbl.', ptID, ch);
            Sbl_all{ch} = [];
        else
            baselineMean = baselineSum ./ baselineCount;
    
            Sbl = S ./ repmat(baselineMean, 1, size(S, 2), size(S, 3));

            % keep from choiceFeedbackWindow to the end
            Sbl = Sbl(:, choiceFeedbackWindow:end, :);

            Sbl_all{ch} = Sbl;
        end

        clear S baselineSum baselineMean W bW pow bPow Sbl

    end

    % Save participant data

    freqSaved = freq_all{1,1};

    saveFileName = fullfile(output_folder, [ptID '_spectrogram_data.mat']);

    save(saveFileName, ...
        'Sbl_all', ...
        'eventTimes', ...
        'bhvData', ...
        'freqSaved', ...
        'selectedChans', ...
        'SelectedAnatomicalLoc', ...
        '-v7.3');

    clear Sbl_all freq_all freqSaved LFPData NS2 bhvData eventTimes

end

