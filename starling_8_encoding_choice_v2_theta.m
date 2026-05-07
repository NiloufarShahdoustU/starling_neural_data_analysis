% this is correct version. and the same v3 version of outcome but this time
% for choice
% freqBand = [80 120]; for checking if they are encoding outcome.
clc;
clear;
close all;

differentPatients = {'202514', '202518', '202521', '202522' , '202601'};

%% all these times are in ms
cueStart = 200;
cueEnd = 800;

flipStart = 500;
flipEnd = 500;

choiceStart = 1500-1;
choiceEnd = 1000;

choiceFeedbackStart = 1000;
choiceFeedbackEnd = 1500;

totalRewardStart = 200;
totalRewardEnd = 800;

cardShowWindow = 1000;
choiceFeedbackWindow = 1500;

%% Main folder
output_folder = fullfile('\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_8_encoding_choice_v2\');
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
%patients

% 
% ptNumber = 1;
% for p = ptNumber:ptNumber
for p = 1:numel(ptIDs)
    ptID = ptIDs{p};
    fprintf('\n--- reading ptID: %s ---\n', ptID);
    
    input_folder_pt = fullfile(input_folder, ptID); 

    output_folder_pt = fullfile(output_folder, ptID); 
    if ~exist(output_folder_pt, 'dir')
        mkdir(output_folder_pt);
    end

    %% reading neural data
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
    
    %% ns2 
    [nevPath,nevName,~] = fileparts(nevFile);
    NS2 = openNSx(fullfile(nevPath,[nevName '.ns2']));
    original_freq = NS2.MetaTags.SamplingFreq;

    %% reading selected channels using ptTrodesStarling that uses Electrodes.mat
    [trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesSTARLING(ptID);
    selectedChans = find(~isECoG);
    selectedChans = selectedChans(1:end-1); 
    SelectedAnatomicalLoc = anatomicalLocs(selectedChans); 

        % --- Remove any anatomical locations that start with "NAC" ---
    maskKeep = ~startsWith(SelectedAnatomicalLoc, 'NAC', 'IgnoreCase', true);
    selectedChans = selectedChans(maskKeep);
    SelectedAnatomicalLoc = SelectedAnatomicalLoc(maskKeep);



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
            clear tmpData tmp
        else
            tmp = filtfilt(b1, a1, double(NS2.Data(selectedChans(ch),:)));
            tmp = filtfilt(b2, a2, tmp);
            LFPData(ch,:) = tmp;
            clear tmp
        end  
    end

    %% now let's denoise data using Elliot's function:
    LFPData = remove1stPC(LFPData);

    %% now let's denoise data using common average rereference:
    % LFPData = LFPData - mean(LFPData, 2);

    %% reading task data:
    bhvFiles = dir(fullfile(input_folder_pt, 'task_data*.csv'));
    if isempty(bhvFiles), error('No task_data*.csv found for %s', ptID); end
    bhvFile = fullfile(bhvFiles(1).folder, bhvFiles(1).name);
    bhvData = readtable(bhvFile);
    
    % remove all rows where trialtype is 'timeout'
    if any(strcmpi(bhvData.Properties.VariableNames, 'trialType'))
        bhvData(strcmpi(bhvData.trialType, 'timeout'), :) = [];
    end
    
    nTrials = size(bhvData, 1);

    %% reading eventTimes:
    eventTimesFile = fullfile(input_folder_pt, 'eventTimes.mat');
    eventTimes = load(eventTimesFile);
    choiceOutcome = eventTimes.choiceAndFeedbackTime;
    trialStart = eventTimes.trialStartTime;

    %% model behavior outputs
    model_root = '\\155.100.91.44\d\Code\Nill\Starling_neural_data\models_outputs';
    model_folders = { ...
        '13_RL_agent_TDlearn_output_greedy', ...
        '13_RL_agent_TDlearn_output_risk_sensitive', ...
        '13_RL_agent_TDlearn_output_softmax'};
    
    model_labels = {'greedy','rs','softmax'};
    model_data = struct();
    
    %% models data handling
    [~, taskName, ~] = fileparts(bhvFile);
    suffix = erase(taskName, 'task_data');

    for m = 1:numel(model_folders)
        model_dir = fullfile(model_root, model_folders{m}, 'model_behavior');
        modelFile = fullfile(model_dir, sprintf('model_behavior%s.csv', suffix));
        if ~exist(modelFile,'file'), error('Model file not found: %s', modelFile); end
    
        % ---- Read as raw text since q_val contains commas ----
        fid = fopen(modelFile,'r');
        header = strsplit(strtrim(fgetl(fid)), ','); 

        textLines = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        lines = textLines{1};
        nT = numel(lines);
    
        % Preallocate
        model_choices = nan(nT,1);
        participant_choices = nan(nT,1);
        model_total_reward = nan(nT,1);
        participant_total_reward = nan(nT,1);
        q_val = cell(nT,1);
        prediction_errors = nan(nT,1);
    
        for t = 1:nT
            L = strtrim(lines{t});
            % Find the first four commas (fields before q_val)
            idx = strfind(L, ',');
            if numel(idx) < 5, continue; end
    
            % Extract fields manually
            parts = cell(1,6);
            parts{1} = strtrim(L(1:idx(1)-1));
            parts{2} = strtrim(L(idx(1)+1:idx(2)-1));
            parts{3} = strtrim(L(idx(2)+1:idx(3)-1));
            parts{4} = strtrim(L(idx(3)+1:idx(4)-1));
    
            % q_val is between 5th comma and last comma
            lastComma = idx(end);
            parts{5} = strtrim(L(idx(4)+1:lastComma-1));
            parts{6} = strtrim(L(lastComma+1:end));
    
            % Parse numerics
            model_choices(t)            = str2double(parts{1});
            participant_choices(t)      = str2double(parts{2});
            model_total_reward(t)       = str2double(parts{3});
            participant_total_reward(t) = str2double(parts{4});
            prediction_errors(t)        = str2double(parts{6});
    
            % Parse q_val [[[...]]] array
            nums = regexp(parts{5}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match');
            vals = str2double(nums);
            if numel(vals) == 9*3*2
                q_val{t} = reshape(vals, [9,3,2]);
            else
                q_val{t} = nan(9,3,2);
            end
        end
    
        % Make table
        rawTbl = table(model_choices, participant_choices, ...
            model_total_reward, participant_total_reward, ...
            q_val, prediction_errors);
    
        % Compute regressors
        nT = height(rawTbl);
        Q_up   = nan(nT,1);
        Q_down = nan(nT,1);
        RPE    = rawTbl.prediction_errors;
    
        deckMap = containers.Map({'uniform','low','high'},{1,2,3});
        cardNum  = bhvData.myCard;
        deckType = string(bhvData.distribution);
    
        for t = 1:nT
            qMat = rawTbl.q_val{t,1}; % 9x3x2
            if any(isnan(qMat),'all') || isnan(cardNum(t)) || ~isKey(deckMap, deckType(t))
                Q_up(t) = nan; Q_down(t) = nan;
                continue
            end
            c = cardNum(t);
            d = deckMap(deckType(t));
            Q_down(t) = qMat(c,d,2);
            Q_up(t)   = qMat(c,d,1);
        end
    
        dQ     = Q_up - Q_down;
        abs_dQ = abs(dQ);
        RPEpos = max(RPE,0);
        RPEneg = min(RPE,0);
    
        modelTbl = rawTbl;
        modelTbl.Q_up   = Q_up;
        modelTbl.Q_down = Q_down;
        modelTbl.dQ     = dQ;
        modelTbl.abs_dQ = abs_dQ;
        modelTbl.RPE    = RPE;
        modelTbl.RPEpos = RPEpos;
        modelTbl.RPEneg = RPEneg;
    
        model_data.(model_labels{m}) = modelTbl;
    end

    %% spectogram:
    fprintf('\nComputing spectrograms and encoding for %s...\n', ptID);
    fWin = [1 200];
    waitBar = 0;
    motherWaveletParam = 6;
    
    S_all = cell(1, nChans);
    baseLineS_all = cell(1, nChans);
    freq_all = cell(1, nChans);
    
    % build once: valid trials with known choiceOutcome
    validTrials_global = find(~isnan(choiceOutcome));

    for ch = 1:nChans
        fprintf('\nDoing spectral calculations for chan %d of %d', ch, nChans)

        % RESET per channel to avoid contamination
        S = [];
        baseLineS = [];

        for ii = 1:numel(validTrials_global)
            tt = validTrials_global(ii);
    
            % LFP segment: choice-centered window [-choiceFeedbackWindow, +choiceFeedbackWindow-1]
            whichData = (choiceOutcome(tt) - choiceFeedbackWindow):(choiceOutcome(tt) + choiceFeedbackWindow- 1);
            whichData = whichData(whichData > 0 & whichData <= size(LFPData,2));
            LFPseg = LFPData(ch, whichData);
    
            [W, period, ~] = basewaveERP(LFPseg, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
            S(:, :, ii) = 10*log10(abs(W).^2);
    
            if ii == 1
                freq_all{ch} = 1 ./ period;
            end
    
            % baseline (-750 to +500 ms around trialStart)
            baseIdx = (trialStart(tt) - 750):(trialStart(tt) + 500 - 1);
            baseIdx = baseIdx(baseIdx > 0 & baseIdx <= size(LFPData,2));
            LFPbase = LFPData(ch, baseIdx);
            [bW, ~, ~] = basewaveERP(LFPbase, original_freq, fWin(1), fWin(2), motherWaveletParam, waitBar);
            baseLineS(:, :, ii) = 10*log10(abs(bW).^2);
        end
    
        S_all{ch} = S;
        baseLineS_all{ch} = baseLineS;
    end
    %spectogram channels end
    
    %% baseline normalize
    Sbl_all = cell(1, nChans);
    for ch = 1:nChans
        Spec = S_all{ch};
        baseLineS = baseLineS_all{ch};
        % baseline normalize:
        baseline_mean = mean(mean(baseLineS(:,100:end-150), 2), 3, 'omitnan');
        Sbl_all{ch} = Spec ./ baseline_mean;  
    end
    %% extract mean power in high-frequency band 
    freqBand = [4 7];
    power_all = nan(nChans, nTrials);
    
    for ch = 1:nChans
        Spec = Sbl_all{ch};
        if isempty(Spec), continue; end

        freqs = freq_all{ch};
        fIdx = freqs >= freqBand(1) & freqs <= freqBand(2);

        % safe time window (index "zero" ~ choiceFeedbackWindow)
        T = size(Spec,2);
        prefWindow = max(1, choiceFeedbackWindow-500):min(T, choiceFeedbackWindow-1);

        % Average over frequency and time only → keep single-trial power values (valid trials only)
        p_valid = squeeze(mean(mean(Spec(fIdx, prefWindow, :), 1, 'omitnan'), 2, 'omitnan')); % [nValidTrials x 1]
        
        % map back to full trial space
        power_all(ch, validTrials_global) = p_valid(:).';
    end

    %% debug
    assert(height(model_data.greedy) == nTrials, 'greedy rows != nTrials');
    assert(height(model_data.softmax) == nTrials, 'softmax rows != nTrials');
    assert(height(model_data.rs)     == nTrials, 'rs rows != nTrials');

    %% find out significant channels
    alpha_level = 0.05;  % significance threshold

    sigChans = false(1, nChans);
    pvals = nan(1, nChans);

    for ch = 1:nChans
        Spec = Sbl_all{ch};
        if isempty(Spec)
            pvals(ch) = NaN; continue
        end

        freqs = freq_all{ch};
        fIdx = freqs >= freqBand(1) & freqs <= freqBand(2);
        bandPower = squeeze(mean(Spec(fIdx, :, :), 1, 'omitnan'));  % time × trials

        % define in-window (prefWindow) and out-of-window indices
        allTimeIdx = 1:size(bandPower, 1);
        T = numel(allTimeIdx);

        outWindowIdx = setdiff(allTimeIdx, prefWindow);

        % average power over trials for each timepoint
        meanPowerTime = mean(bandPower, 2, 'omitnan');

        % sample in-window vs. out-window
        inPower  = meanPowerTime(prefWindow);
        outPower = meanPowerTime(outWindowIdx);

        if isempty(inPower) || isempty(outPower)
            pvals(ch) = NaN;
        else
            [~, p] = ttest2(inPower, outPower);
            pvals(ch) = p;
        end
    end

    sigIdx = find(pvals < alpha_level);
    if isempty(sigIdx)
        warning('No channels passed HFB screen; proceeding with all ECoG channels.');
        sigIdx = 1:nChans;
    end

    
    NewAnatomicalLocs = SelectedAnatomicalLoc(sigIdx);
    
    %% ridge-regression encoding 
    model_list = fieldnames(model_data);
    lambda     = 1;      % fixed ridge penalty
    minN       = 40;     % skip channels with too few valid trials
    alpha_level = 0.05;  % significance threshold
    
    results = struct();
    
    for m = 1:numel(model_list)
        mdl = model_list{m};
        fprintf('Encoding for %s model...\n', mdl);
    
        % ---------- Build regressors ----------
        if strcmp(mdl, 'rs')
            Xraw = [model_data.rs.dQ, ...
                    model_data.rs.abs_dQ, ...
                    model_data.rs.RPE, ...
                    model_data.rs.RPEpos, ...
                    model_data.rs.RPEneg];
        else
            Xraw = [model_data.(mdl).dQ, ...
                    model_data.(mdl).abs_dQ, ...
                    model_data.(mdl).RPE];
        end
    
        % ---------- Preallocate ----------
        cvPseudoR2 = nan(nChans,1);
        pvals      = nan(nChans,1);
    
        % ---------- Loop over channels ----------
        for ch = sigIdx
            % --- Response variable ---
            yfull = power_all(ch,:)';   % [nTrials x 1]
    
            % --- Align lengths ---
            nXY  = min(size(Xraw,1), numel(yfull));
            Xuse = Xraw(1:nXY, :);
            yuse = yfull(1:nXY);
    
            % --- Remove NaNs ---
            valid = all(isfinite(Xuse),2) & isfinite(yuse);
            Xv = Xuse(valid,:);
            yv = yuse(valid);
    
            if numel(yv) < minN
                continue
            end
    
            % --- Add intercept ---
            Xv = [ones(size(Xv,1),1), Xv];
            I = eye(size(Xv,2)); I(1,1) = 0; % no penalty on intercept
    
            % --- Ridge regression ---
            beta = (Xv'*Xv + lambda*I) \ (Xv'*yv);
            yhat = Xv * beta;
    
            % --- Compute pseudo-R² ---
            SSE = nansum((yv - yhat).^2);
            SST = nansum((yv - mean(yv)).^2);
            cvPseudoR2(ch) = 1 - SSE/max(SST, eps);
    
            % --- Analytic F-test for significance ---
            n = numel(yv);
            k = size(Xv,2) - 1; % exclude intercept
            R2 = max(0, cvPseudoR2(ch));
            Fval = (R2 / k) / ((1 - R2) / (n - k - 1));
            pvals(ch) = 1 - fcdf(Fval, k, n - k - 1);
        end
    
        % --- Significance mask ---
        sigMask = pvals < alpha_level;
    
        % --- Store results ---
        results.(mdl).cvPseudoR2 = cvPseudoR2;
        results.(mdl).pvals      = pvals;
        results.(mdl).sigMask    = sigMask;
        results.(mdl).alpha_level = alpha_level;
        results.(mdl).lambda     = lambda * ones(size(cvPseudoR2));
    end
    
    
    model_list = fieldnames(results);
    nChans_total = [];
    keepMask = [];
    
    for m = 1:numel(model_list)
        if isfield(results.(model_list{m}), 'cvPseudoR2')
            v = results.(model_list{m}).cvPseudoR2(:);
            nChans_total = numel(v);
            keepMask = false(nChans_total,1);
            break
        end
    end
    
    if ~isempty(nChans_total)
        for m = 1:numel(model_list)
            if isfield(results.(model_list{m}), 'cvPseudoR2')
                v = results.(model_list{m}).cvPseudoR2(:);
                if numel(v) == nChans_total
                    keepMask = keepMask | isfinite(v);
                end
            end
        end
    
        if exist('pvals','var') && numel(pvals) == nChans_total
            keepMask = keepMask & isfinite(pvals(:));
        end
    
        for m = 1:numel(model_list)
            fns = fieldnames(results.(model_list{m}));
            for fi = 1:numel(fns)
                vec = results.(model_list{m}).(fns{fi});
                if isnumeric(vec) && isvector(vec) && numel(vec) == nChans_total
                    results.(model_list{m}).(fns{fi}) = vec(keepMask);
                elseif iscell(vec) && isvector(vec) && numel(vec) == nChans_total
                    results.(model_list{m}).(fns{fi}) = vec(keepMask);
                end
            end
        end
    
        if exist('NewAnatomicalLocs','var') && numel(NewAnatomicalLocs) == nChans_total
            if iscell(NewAnatomicalLocs) || isstruct(NewAnatomicalLocs)
                NewAnatomicalLocs = NewAnatomicalLocs(keepMask);
            end
        end
    
        if exist('pvals','var') && numel(pvals) == nChans_total
            pvals = pvals(keepMask);
        end
    
        if exist('sigIdx','var') && ~isempty(sigIdx)
            old2new = zeros(nChans_total,1);
            old2new(keepMask) = 1:nnz(keepMask);
            sigIdx = old2new(sigIdx(:));
            sigIdx = sigIdx(sigIdx > 0);
        end
    end
    
    save(fullfile(output_folder_pt, 'encodingResults.mat'), ...
         'results','NewAnatomicalLocs','freqBand','sigIdx','alpha_level');


%% vis

        
    labels = {'greedy', 'softmax', 'RS'};
    hexColors = {'#56B4E9', '#009E73', '#065B8D'};  % your color palette
    colors = cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x',[1 3])/255, hexColors, 'UniformOutput', false);
    
    % ====== Only keep significant electrodes (e.g., p < 0.05) ======
    sigThresh = 0.05;
    
    data = cell(1,3);
    data{1} = results.greedy.cvPseudoR2(results.greedy.pvals < sigThresh);
    data{2} = results.softmax.cvPseudoR2(results.softmax.pvals < sigThresh);
    data{3} = results.rs.cvPseudoR2(results.rs.pvals < sigThresh);
    
    for i = 1:3
        data{i} = data{i}(~isnan(data{i}));
    end
    
    % --- Combine data into one vector and a grouping variable ---
    allData = vertcat(data{:});
    group = [repmat({'greedy'}, numel(data{1}), 1);
             repmat({'softmax'}, numel(data{2}), 1);
             repmat({'RS'}, numel(data{3}), 1)];
    
    figure('Color','w','Position',[100 100 600 400]); hold on;
    
    % --- Boxplot (grouped format) ---
    boxplot(allData, group, ...
        'Colors', 'k', ...
        'Symbol', '', ...
        'BoxStyle', 'outline');
    
    % --- Overlay scatter points per group ---
    for i = 1:3
        x = i * ones(size(data{i}));
        scatter(x + randn(size(x))*0.05, data{i}, 20, colors{i}, 'filled', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
    end
    
    ylabel('R^2 per electrode (significant only)');
    title(sprintf('Significant encoding power – %s', ptID), 'FontWeight','bold');
    set(gca, 'Box', 'off', 'FontSize', 12);
    axis square;
    
    % ===== Dynamic y-limits =====
    if isempty(allData)
        allData = 0;
    end
    dmin = min(allData);
    dmax = max(allData);
    drng = dmax - dmin;
    if drng == 0, drng = max(0.1, abs(dmax)*0.1 + 1e-6); end
    ylim([dmin - 0.05*drng, dmax + 0.1*drng]);
    
    saveas(gcf, fullfile(output_folder_pt, [ptID '_encoding_significant.pdf']));
    close(gcf);


end



%% debug

% assert(height(model_data.greedy) == nTrials, 'greedy rows != nTrials');
% assert(height(model_data.softmax) == nTrials, 'softmax rows != nTrials');
% assert(height(model_data.rs)     == nTrials, 'rs rows != nTrials');
