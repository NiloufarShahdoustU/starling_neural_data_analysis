% outcome based spectogram

% trialStart: 1 , by fixation point: 500ms
% cardShow: 2 , 1000ms
% instructionMessage: 3 , no delay
% flipSpace: 4 , no delay <3000ms
% choiceAndFeedback: 5 , feedback is shown for 2000ms
% totalReward: 6 , shown for 1000ms

clc;
clear;
close all;
differentPatients = {'202421', '202509', '202511', '202512'};
needsConcatenationPts = {'202514','202521'};
%% all these times are in ms
cueStart = 200;
cueEnd = 800;

flipStart = 500;
flipEnd = 500;

choiceFeedbackStart = 1000;
choiceFeedbackEnd = 1500;

totalRewardStart = 200;
totalRewardEnd = 800;

cardShowWindow = 1000;
choiceFeedbackWindow = 1500;

%% Main folder
output_folder = fullfile('\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_7_spectogram_outcome\');
input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');
d = dir(input_folder);

isub = [d(:).isdir]; 
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = string(subFolders);

%%

% ptNumber = 1;
%patients
% for p = ptNumber:ptNumber
for p = 1:numel(ptIDs)
    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);
    
    input_folder_pt = fullfile(input_folder, ptID); 
    output_folder_pt = fullfile(output_folder, ptID); 

   % reading neural data
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
        if ismember(ptID, needsConcatenationPts)
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
    % freq_all = cell(1, nChans);
    % 
    % for ch = 1:nChans
    %     fprintf('\nDoing spectral calculations for chan %d of %d', ch, nChans)
    %     for tt = 1:nTrials
    %         % LFP 
    %         if ~isnan(choiceOutcome(tt))
    %             whichData = (choiceOutcome(tt) - choiceFeedbackWindow):(choiceOutcome(tt) + choiceFeedbackWindow- 1);
    %             LFPseg = LFPData(ch, whichData);
    % 
    %             [W, period, scale] = basewaveERP(LFPseg, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
    %             S(:, :, tt) = 10*log10(abs(W).^2);
    % 
    %             if tt == 1
    %                 freq_all{ch} = 1 ./ period; 
    %             end
    % 
    %             % baseline (-700 to +300 around trialStart)
    %             baseData = (trialStart(tt) - 750):(trialStart(tt) + 500 - 1);
    %             LFPbase = LFPData(ch, baseData);
    %             [bW, bperiod, bscale] = basewaveERP(LFPbase, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
    %             baseLineS(:, :, tt) = 10*log10(abs(bW).^2);
    %         end
    %     end
    % 
    %     % save results for channel
    %     S_all{ch} = S;
    %     baseLineS_all{ch} = baseLineS;
    % end


    freq_all = cell(1, nChans);
    
    for ch = 1:nChans
        fprintf('\nDoing spectral calculations for chan %d of %d', ch, nChans)
    
        clear S
        S = [];
        baselineSum = [];
        baselineCount = 0;
    
        for tt = 1:nTrials
    
            if ~isnan(choiceOutcome(tt))
    
                % main LFP segment
                whichData = (choiceOutcome(tt) - choiceFeedbackWindow):(choiceOutcome(tt) + choiceFeedbackWindow- 1);
                LFPseg = LFPData(ch, whichData);
    
                [W, period, scale] = basewaveERP(LFPseg, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
                pow = single(10*log10(abs(W).^2));
    
                S(:, :, tt) = pow;
    
                if isempty(freq_all{ch})
                    freq_all{ch} = 1 ./ period;
                end
    
                % baseline segment
                baseData = (trialStart(tt) - 750):(trialStart(tt) + 500 - 1);
                LFPbase = LFPData(ch, baseData);
    
                [bW, bperiod, bscale] = basewaveERP(LFPbase, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
                bPow = single(10*log10(abs(bW).^2));
    
                baselineWindow = bPow(:, 100:end-150);
    
                if isempty(baselineSum)
                    baselineSum = zeros(size(bPow,1), 1, 'single');
                end
    
                baselineSum = baselineSum + mean(baselineWindow, 2);
                baselineCount = baselineCount + 1;
            end
        end
    
        baselineMean = baselineSum ./ baselineCount;
    
        Sbl_all{ch} = S ./ repmat(baselineMean, 1, size(S,2), size(S,3));
    
        clear S baselineSum baselineMean W bW pow bPow
    end


    
    % % baseline normalize 
    % Sbl_all = cell(1, nChans);
    % for ch = 1:nChans
    %     Spec = S_all{ch};
    %     baseLineS = baseLineS_all{ch};
    %     % baseline normalize:
    %     Sbl_all{ch} = Spec ./ repmat(mean(mean(baseLineS(:,100:end-150), 2), 3), 1, size(S, 2), size(S, 3));
    % 
    %     % frequency noramlize:
    %     % Sbl_all{ch} = S ./ repmat(period', 1, size(S, 2), size(S, 3));
    % end
    
    
    winTrials = find(bhvData.outcome == "win");
    loseTrials = find(bhvData.outcome == "lose");


    
    condNames = {'win', 'lose'};
    condTrials = {winTrials, loseTrials};

    % % Count number of win/lose trials
    % nWin = numel(winTrials);
    % nLose = numel(loseTrials);
    % nTotal = nWin + nLose;
    % 
    % % Normalization factors (so that both are scaled to same total)
    % winNormFactor = nWin / nTotal;
    % loseNormFactor = nLose / nTotal;

    
    % visualization
    fig_folder = fullfile(output_folder_pt, 'spectrogram');
    if ~exist(fig_folder, 'dir')
        mkdir(fig_folder);
    end
    
    nRows = 10;
    nCols = ceil(nChans / nRows);
    
    for c = 1:numel(condNames)
        trialsIdx = condTrials{c};
    
        % colorbar scale across all channels
        allVals = [];
        for ch = 1:nChans
            tmp = mean(Sbl_all{ch}(:, choiceFeedbackWindow-choiceFeedbackStart:choiceFeedbackWindow+choiceFeedbackEnd-1, trialsIdx), 3);
            allVals = [allVals; tmp(:)];
        end
        climVals = [prctile(allVals, 5), prctile(allVals, 95)];
    
        
        f = figure('Visible', 'off', 'Position', [100 100 1000 1000]);
        sgtitle(sprintf('%s | %s', ptID, condNames{c}), 'FontWeight', 'bold');
    
        for ch = 1:nChans
            subplot(nRows, nCols, ch);
    
            dataToPlot = mean(Sbl_all{ch}(:, choiceFeedbackWindow-choiceFeedbackStart:choiceFeedbackWindow+choiceFeedbackEnd-1, trialsIdx), 3);

            % % Normalize by condition trial count proportion
            % if strcmp(condNames{c}, 'win')
            %     dataToPlot = dataToPlot ./ winNormFactor;
            % else
            %     dataToPlot = dataToPlot ./ loseNormFactor;
            % end


            imagesc(1:size(dataToPlot,2), freq_all{ch}, dataToPlot);
            set(gca, 'YDir', 'normal', 'YScale', 'log');

            ylim([fWin(1) fWin(2)]);
            % ----------------------------------------------------------
            
            axis square tight;
            caxis(climVals); % same colorbar scale for all
            set(gca, 'FontSize', 2);
    
            %  red vertical line at t = 200 ms where back of the cards appear
            hold on;
            xline(1000, 'r', 'LineWidth', 0.5);
            hold off;


            title(sprintf('%s', anatomicalLocs{selectedChans(ch)}), ...
      'FontSize', 3, 'FontWeight', 'normal', 'Interpreter', 'none');
        end
    
        %  one shared horizontal colorbar at the bottom
        h = colorbar('southoutside');
        h.Position = [0.25 0.05 0.5 0.02];
        h.Label.String = 'power (normalized)';
        h.FontSize = 8;
        
        set(f,'Renderer','painters');
        exportgraphics(f, ...
            fullfile(fig_folder, sprintf('%s_%s.pdf', ptID, condNames{c})), ...
            'ContentType','vector', ...
            'BackgroundColor','none', ...
            'Resolution',600);
        close(f);
    end


end
%patients

%% DEEEEEEBBBBUUUUUGGGG


% close all
% dataToPlot = mean(baseLineS_all{1}(:, :, 1), 3);
% imagesc(dataToPlot);
% set(gca, 'YDir', 'normal');

% 
% aaa = 1:size(dataToPlot,2);
% close all
% plot(log10(freq_all{1}));


