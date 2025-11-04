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
differentPatients = {'202514', '202518'};

%% all these times are in ms
cueStart = 200;
cueEnd = 800;

flipStart = 1500-1;
flipEnd = 500;

choiceFeedbackStart = 1000;
choiceFeedbackEnd = 1500;

totalRewardStart = 200;
totalRewardEnd = 800;

cardShowWindow = 1000;

flipWindow = 1500;

%% Main folder
output_folder = fullfile('\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_7_spectogram_flip\');
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
    selectedChans = find(isECoG);
    selectedChans = selectedChans(1:end-1); 
    SelectedAnatomicalLoc = anatomicalLocs(selectedChans);

    nChans = length(selectedChans);
    LFPData = [];
    [b1,a1] = iirnotch(60/(original_freq/2), (60/(original_freq/2))/25);
    [b2,a2] = iirnotch(120/(original_freq/2), (120/(original_freq/2))/25);
    
    for ch = 1:nChans
        if strcmp(ptID, '202514')
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
    flipTimes = eventTimes.flipSpaceTime;
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
            if ~isnan(flipTimes(tt))
                whichData = (flipTimes(tt) - flipWindow):(flipTimes(tt) + flipWindow- 1);
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

        % frequency noramlize:
        % Sbl_all{ch} = S ./ repmat(period', 1, size(S, 2), size(S, 3));
    end
    
    
    UniTrials = find(bhvData.distribution == "uniform");
    LowTrials = find(bhvData.distribution == "low");
    HighTrials = find(bhvData.distribution == "high");

    
    condNames = {'Uniform', 'Low', 'High'};
    condTrials = {UniTrials, LowTrials, HighTrials};


    
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
            tmp = mean(Sbl_all{ch}(:, flipWindow-flipStart:flipWindow+flipEnd-1, trialsIdx), 3);
            allVals = [allVals; tmp(:)];
        end
        climVals = [prctile(allVals, 5), prctile(allVals, 95)];
    
        
        f = figure('Visible', 'off', 'Position', [100 100 1000 1000]);
        sgtitle(sprintf('%s | %s', ptID, condNames{c}), 'FontWeight', 'bold');
    
        for ch = 1:nChans
            subplot(nRows, nCols, ch);
    
            dataToPlot = mean(Sbl_all{ch}(:, flipWindow-flipStart+1:flipWindow+flipEnd-1, trialsIdx), 3);


            imagesc(1:size(dataToPlot,2), freq_all{ch}, dataToPlot);
            set(gca, 'YDir', 'normal', 'YScale', 'log');

            ylim([fWin(1) fWin(2)]);
            % ----------------------------------------------------------
            
            axis square tight;
            caxis(climVals); % same colorbar scale for all
            set(gca, 'FontSize', 2);
    
            %  red vertical line at t = 200 ms where back of the cards appear
            hold on;
            xline(1500, 'r', 'LineWidth', 0.5);
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


