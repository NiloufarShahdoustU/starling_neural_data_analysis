% outcome based spectogram saved!!









% this code might not be correct!!!!!
clc;
clear;
close all;
close all;
differentPatients = {'202514', '202518', '202521', '202522' , '202601'};

%% all these times are in ms
cueStart = 200;
cueEnd = 800;

flipStart = 500;
flipEnd = 500;

choiceStart = 1000;
choiceEnd = 1000;

choiceFeedbackStart = 1000;
choiceFeedbackEnd = 1500;

totalRewardStart = 200;
totalRewardEnd = 800;

cardShowWindow = 1000;
choiceFeedbackWindow = 1500;

%% Main folder
output_folder = fullfile('\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_9_spectogram_4_areas_prep\');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');
d = dir(input_folder);

isub = [d(:).isdir]; 
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = string(subFolders);

%%

ptNumber = 9;
for p = ptNumber:ptNumber
% for p = 1:numel(ptIDs)
    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);
    
        
    input_folder_pt = fullfile(input_folder, ptID);  


   % reading neural data
    if any(strcmp(ptID, differentPatients ))
        nevList = dir(fullfile(input_folder_pt, '\hub_neural_data\*.nev'));
    else
        nevList = dir(fullfile(input_folder_pt, '*.nev'));
    end
    
    if length(nevList) > 1
        error('many nev files available for this patient. Please specify...')
    elseif isempty(nevList)
        error('no nev files found...')
    else
        nevFile = fullfile(nevList.folder, nevList.name);
    end
    
    % ns2 
    [nevPath,nevName,~] = fileparts(nevFile);
    NS2 = openNSx(fullfile(nevPath,[nevName '.ns2']));
    original_freq = NS2.MetaTags.SamplingFreq;

    % reading selected channels using ptTrodesStarling that uses
    % Electrodes.mat
    [trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesSTARLING(ptID);
    selectedChans = find(~isECoG);
    selectedChans = selectedChans(1:end-1); 
    SelectedAnatomicalLoc = anatomicalLocs(selectedChans);

    nChans = length(selectedChans);
    LFPData = [];
    [b1,a1] = iirnotch(60/(original_freq/2), (60/(original_freq/2))/25);
    [b2,a2] = iirnotch(120/(original_freq/2), (120/(original_freq/2))/25);
    
    for ch = 1:nChans
        if ismember(ptID, {'202514','202521'})
            % concatenate all cells in NS2.Data along time dimension
            tmpData = [];
            for c = 1:numel(NS2.Data)
                tmpData = [tmpData, double(NS2.Data{c}(selectedChans(ch),:))];
            end
    
            % apply both notch filters sequentially
            tmp = filtfilt(b1, a1, tmpData);
            tmp = filtfilt(b2, a2, tmp);
            LFPData(ch,:) = tmp;
            clear tmpData
    
        else
            tmp = filtfilt(b1, a1, double(NS2.Data(selectedChans(ch),:)));
            tmp = filtfilt(b2, a2, tmp);
            LFPData(ch,:) = tmp;
        end  
        clear tmp
    end

    % now let's denoise data using Elliot's function:
    LFPData = remove1stPC(LFPData);

    % now let's denoise data using common average rereference:
    % LFPData = LFPData - mean(LFPData, 2);

   

   % reading task data:
    bhvFiles = dir(fullfile(input_folder_pt, 'task_data*.csv'));
    bhvFile = fullfile(bhvFiles(1).folder, bhvFiles(1).name);
    bhvData = readtable(bhvFile);
    nTrials = size(bhvData, 1);


    % reading eventTimes:
    eventTimesFile = fullfile(input_folder_pt, 'eventTimes.mat');
    eventTimes = load(eventTimesFile);
    choiceOutcome = eventTimes.choiceAndFeedbackTime;
    trialStart = eventTimes.trialStartTime;
    
    % spectogram

    fWin = [1 200];
    waitBar = 0;
    motherWaveletParam = 6;

    %channels
    S_all = cell(1, nChans);
    baseLineS_all = cell(1, nChans);
    freq_all = cell(1, nChans);

    for ch = 1:nChans
        fprintf('\nDoing spectral calculations for chan %d of %d', ch, nChans)
        for tt = 1:nTrials
            % LFP 
            if ~isnan(choiceOutcome(tt))
                whichData = (choiceOutcome(tt) - choiceFeedbackWindow):(choiceOutcome(tt) + choiceFeedbackWindow- 1);
                LFPseg = LFPData(ch, whichData);
        
                [W, period, scale] = basewaveERP(LFPseg, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
                S(:, :, tt) = 10*log10(abs(W).^2);
                
                if tt == 1
                    freq_all{ch} = 1 ./ period; 
                end
        
                % baseline (-700 to +300 around trialStart)
                baseData = (trialStart(tt) - 750):(trialStart(tt) + 500 - 1);
                LFPbase = LFPData(ch, baseData);
                [bW, bperiod, bscale] = basewaveERP(LFPbase, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
                baseLineS(:, :, tt) = 10*log10(abs(bW).^2);
            end
        end
    
        % save results for channel
        S_all{ch} = S;
        baseLineS_all{ch} = baseLineS;
    end
    
    % baseline normalize 
    Sbl_all = cell(1, nChans);
    for ch = 1:nChans
        Spec = S_all{ch};
        baseLineS = baseLineS_all{ch};
        % baseline normalize:
        Sbl_all{ch} = Spec ./ repmat(mean(mean(baseLineS(:,100:end-150), 2), 3), 1, size(S, 2), size(S, 3));
    end
    

%%

    InterestTrials = find(bhvData.choice == "arrowup" | bhvData.choice == "arrowdown");


    dataToSave = nan(nChans,length(period), 2*choiceStart);
    freqToSave = freq_all{1,1};

    for ch = 1:nChans
         dataToSave(ch,:,:) = mean(Sbl_all{ch}(:, choiceFeedbackWindow-choiceStart:choiceFeedbackWindow+choiceEnd-1, InterestTrials), 3);
    end

    % ---------------- Save results for this patient ----------------
    save_path = fullfile(output_folder, ['spectogram_outcome_' ptID '.mat']);
    save(save_path, 'dataToSave', 'SelectedAnatomicalLoc', 'freqToSave', '-v7.3');
    fprintf('\nSaved spectrogram data for %s to:\n%s\n', ptID, save_path);

end
%patients


