clc;
clear;
close all;

output_folder = fullfile('\\\\155.100.91.44\\d\\Code\\Nill\\Starling_neural_data\\starling_12_CBPT_choice_v2\\');
input_folder  = fullfile('\\\\155.100.91.44\\d\\Data\\Nill\\starling\\spectrograms\\choice\\');

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

files = dir(fullfile(input_folder, '*_spectrogram_data.mat'));

nPerm = 1000;

alphaBin = 0.10;
alphaCluster = 0.05;
minClusterSize = 5;
minClusterDurationMs = 50;
minSigPixelRunMs = 40;

fWin = [5 200];
removeLastMs = 100;

applyTimeSmoothing = true;
timeSmoothSigma = 3;

applyFreqSmoothing = true;
freqSmoothSigma = 1;

channelCombineMethod = 'mean';

rng(1);

for p = 1:numel(files)

    file_name = files(p).name;
    file_path = fullfile(files(p).folder, file_name);

    ptID = erase(file_name, '_spectrogram_data.mat');

    fprintf('\nProcessing ptID: %s\n', ptID);

    S = load(file_path);

    Sbl_all        = S.Sbl_all;
    bhvData        = S.bhvData;
    freq_all       = S.freqSaved;
    anatomicalLocs = S.SelectedAnatomicalLoc;

    nChans = numel(Sbl_all);

    choice = lower(string(bhvData.choice));

    arrowUpTrials   = find(choice == "arrowup");
    arrowDownTrials = find(choice == "arrowdown");

    condNames  = {'arrowup', 'arrowdown'};
    condTrials = {arrowUpTrials, arrowDownTrials};

    cleanAnat = strings(nChans, 1);

    for ch = 1:nChans

        thisLoc = string(anatomicalLocs{ch});

        thisLoc = regexprep(thisLoc, 'left', '', 'ignorecase');
        thisLoc = regexprep(thisLoc, 'right', '', 'ignorecase');
        thisLoc = regexprep(thisLoc, '\\blh\\b', '', 'ignorecase');
        thisLoc = regexprep(thisLoc, '\\brh\\b', '', 'ignorecase');
        thisLoc = regexprep(thisLoc, '(^|\\W)L(\\W|$)', ' ', 'ignorecase');
        thisLoc = regexprep(thisLoc, '(^|\\W)R(\\W|$)', ' ', 'ignorecase');

        thisLoc = regexprep(thisLoc, '[-_]', ' ');
        thisLoc = strtrim(regexprep(thisLoc, '\\s+', ' '));

        cleanAnat(ch) = thisLoc;

    end

    brainAreas = unique(cleanAnat, 'stable');
    brainAreas(brainAreas == "") = [];

    for a = 1:numel(brainAreas)

        areaName = brainAreas(a);
        areaChans = find(cleanAnat == areaName);

        safeAreaName = matlab.lang.makeValidName(char(areaName));

        fprintf('\nRunning sign-separated t-test CBPT for area: %s | channels: %d\n', ...
            areaName, numel(areaChans));

        firstCh = areaChans(1);

        if iscell(freq_all)
            freqVec = freq_all{firstCh};
        else
            freqVec = freq_all;
        end

        nFreq = size(Sbl_all{firstCh}, 1);
        nTime = size(Sbl_all{firstCh}, 2);

        validTimeIdx = 1:(nTime - removeLastMs);
        validFreqIdx = find(freqVec >= fWin(1) & freqVec <= fWin(2));

        freqVecValid = freqVec(validFreqIdx);

        statsTrials = [arrowUpTrials; arrowDownTrials];
        statsChoice = choice(statsTrials);

        areaTrialData = nan(numel(validFreqIdx), numel(validTimeIdx), numel(statsTrials));

        for tr = 1:numel(statsTrials)

            trialIdx = statsTrials(tr);

            chanData = nan(numel(validFreqIdx), numel(validTimeIdx), numel(areaChans));

            for k = 1:numel(areaChans)

                ch = areaChans(k);

                tmp = Sbl_all{ch}(:, :, trialIdx);
                tmp = tmp(validFreqIdx, validTimeIdx);

                if applyTimeSmoothing || applyFreqSmoothing
                    tmp = smooth_time_freq_matrix(tmp, timeSmoothSigma, freqSmoothSigma, applyTimeSmoothing, applyFreqSmoothing);
                end

                chanData(:, :, k) = tmp;

            end

            switch lower(channelCombineMethod)

                case 'median'
                    areaTrialData(:, :, tr) = median(chanData, 3, 'omitnan');

                case 'mean'
                    areaTrialData(:, :, tr) = mean(chanData, 3, 'omitnan');

                otherwise
                    error('Unknown channelCombineMethod.');

            end

        end

        arrowUpLocalIdx   = find(statsChoice == "arrowup");
        arrowDownLocalIdx = find(statsChoice == "arrowdown");

        if isempty(arrowUpLocalIdx) || isempty(arrowDownLocalIdx)

            fprintf('Skipping %s because one condition is empty.\n', areaName);

            sigMaskSmall = false(numel(validFreqIdx), numel(validTimeIdx));
            sigMaskPosSmall = false(numel(validFreqIdx), numel(validTimeIdx));
            sigMaskNegSmall = false(numel(validFreqIdx), numel(validTimeIdx));

            pMapSmall = nan(numel(validFreqIdx), numel(validTimeIdx));
            tMapSmall = nan(numel(validFreqIdx), numel(validTimeIdx));

            clusterThresholdPos = NaN;
            clusterThresholdNeg = NaN;
            realClusterTable = table();
            maxNullPos = [];
            maxNullNeg = [];

        else

            [sigMaskSmall, ...
             sigMaskPosSmall, ...
             sigMaskNegSmall, ...
             pMapSmall, ...
             tMapSmall, ...
             clusterThresholdPos, ...
             clusterThresholdNeg, ...
             realClusterTable, ...
             maxNullPos, ...
             maxNullNeg] = run_ttest_cbpt_sign_separated( ...
                areaTrialData, ...
                arrowUpLocalIdx, ...
                arrowDownLocalIdx, ...
                freqVecValid, ...
                validTimeIdx, ...
                nPerm, ...
                alphaBin, ...
                alphaCluster, ...
                minClusterSize, ...
                minClusterDurationMs);

            sigMaskPosSmall = remove_short_time_runs(sigMaskPosSmall, validTimeIdx, minSigPixelRunMs);
            sigMaskNegSmall = remove_short_time_runs(sigMaskNegSmall, validTimeIdx, minSigPixelRunMs);
            sigMaskSmall = sigMaskPosSmall | sigMaskNegSmall;

        end

        fullSigMask = false(nFreq, nTime);
        fullSigMaskPos = false(nFreq, nTime);
        fullSigMaskNeg = false(nFreq, nTime);

        fullPMap = nan(nFreq, nTime);
        fullTMap = nan(nFreq, nTime);

        fullSigMask(validFreqIdx, validTimeIdx) = sigMaskSmall;
        fullSigMaskPos(validFreqIdx, validTimeIdx) = sigMaskPosSmall;
        fullSigMaskNeg(validFreqIdx, validTimeIdx) = sigMaskNegSmall;

        fullPMap(validFreqIdx, validTimeIdx) = pMapSmall;
        fullTMap(validFreqIdx, validTimeIdx) = tMapSmall;

        matName = fullfile(output_folder, ...
            sprintf('%s_%s_choice_ttest_CBPT_results.mat', ...
            ptID, safeAreaName));

        save(matName, ...
            'sigMaskSmall', ...
            'sigMaskPosSmall', ...
            'sigMaskNegSmall', ...
            'fullSigMask', ...
            'fullSigMaskPos', ...
            'fullSigMaskNeg', ...
            'pMapSmall', ...
            'tMapSmall', ...
            'fullPMap', ...
            'fullTMap', ...
            'freqVec', ...
            'freqVecValid', ...
            'validFreqIdx', ...
            'validTimeIdx', ...
            'areaName', ...
            'areaChans', ...
            'condNames', ...
            'nPerm', ...
            'alphaBin', ...
            'alphaCluster', ...
            'minClusterSize', ...
            'minClusterDurationMs', ...
            'minSigPixelRunMs', ...
            'clusterThresholdPos', ...
            'clusterThresholdNeg', ...
            'realClusterTable', ...
            'maxNullPos', ...
            'maxNullNeg', ...
            'applyTimeSmoothing', ...
            'timeSmoothSigma', ...
            'applyFreqSmoothing', ...
            'freqSmoothSigma', ...
            'channelCombineMethod', ...
            '-v7.3');

        if ~isempty(realClusterTable)

            csvName = fullfile(output_folder, ...
                sprintf('%s_%s_significant_cluster_table.csv', ...
                ptID, safeAreaName));

            writetable(realClusterTable, csvName);

        end

        areaCondDataSmoothed = cell(1, 2);
        rowVals = [];

        for c = 1:2

            trialsIdx = condTrials{c};

            if isempty(trialsIdx)
                areaCondDataSmoothed{c} = [];
                continue;
            end

            chanData = [];

            for k = 1:numel(areaChans)

                ch = areaChans(k);

                tmp = mean(Sbl_all{ch}(:, :, trialsIdx), 3, 'omitnan');

                if applyTimeSmoothing || applyFreqSmoothing
                    tmp = smooth_time_freq_matrix(tmp, timeSmoothSigma, freqSmoothSigma, applyTimeSmoothing, applyFreqSmoothing);
                end

                chanData(:, :, k) = tmp;

            end

            switch lower(channelCombineMethod)

                case 'median'
                    areaMeanSmoothed = median(chanData, 3, 'omitnan');

                case 'mean'
                    areaMeanSmoothed = mean(chanData, 3, 'omitnan');

            end

            areaCondDataSmoothed{c} = areaMeanSmoothed;

            rowVals = [rowVals; areaMeanSmoothed(:)];

        end

        rowVals = rowVals(~isnan(rowVals));

        if isempty(rowVals)
            rowClimVals = [];
        else
            rowClimVals = [prctile(rowVals, 5), prctile(rowVals, 95)];
        end

        f = figure('Visible','off', ...
            'Units','inches', ...
            'Position',[1 1 18 9], ...
            'PaperUnits','inches', ...
            'PaperSize',[18 9], ...
            'PaperPosition',[0 0 18 9]);

        sgtitle(sprintf('%s | %s | time-frequency choice spectrogram CBPT', ...
            ptID, areaName), ...
            'FontWeight','bold', ...
            'Interpreter','none');

        axPos = [
            0.07 0.16 0.36 0.70
            0.57 0.16 0.36 0.70
        ];

        for c = 1:2

            ax = axes('Parent', f, ...
                'Position', axPos(c,:));

            dataToPlot = areaCondDataSmoothed{c};

            if isempty(dataToPlot)

                axis off;

                title(sprintf('%s: no trials', condNames{c}), ...
                    'Interpreter','none');

                continue;

            end

            freqMask = freqVec >= fWin(1) & freqVec <= fWin(2);

            timeVec = 1:size(dataToPlot, 2);
            timeVec = timeVec(1:end-removeLastMs);

            dataPlot = dataToPlot(freqMask, 1:end-removeLastMs);
            freqPlot = freqVec(freqMask);

            sigMaskPlot = fullSigMask(freqMask, 1:end-removeLastMs);

            alphaMask = 0.20 * ones(size(sigMaskPlot));
            alphaMask(sigMaskPlot) = 1.0;

            hImg = imagesc(ax, timeVec, freqPlot, dataPlot);

            set(hImg, 'AlphaData', alphaMask);

            set(ax, ...
                'YDir','normal', ...
                'YScale','log', ...
                'FontSize',11, ...
                'LooseInset', ...
                max(get(ax,'TightInset'), 0.04));

            xlim(ax, [timeVec(1), timeVec(end)]);
            ylim(ax, [fWin(1), fWin(2)]);

            if ~isempty(rowClimVals)
                caxis(ax, rowClimVals);
            end

            hold(ax, 'on');
            xline(ax, 1000, 'r', 'LineWidth', 0.5);
            hold(ax, 'off');

            title(ax, sprintf('%s', condNames{c}), ...
                'FontWeight','bold', ...
                'FontSize',14);

            xlabel(ax, 'Time', 'FontSize',12);
            ylabel(ax, 'Frequency', 'FontSize',12);

            cb = colorbar(ax, 'eastoutside');
            cb.Label.String = 'power';
            cb.FontSize = 10;

            box(ax, 'on');

        end

        drawnow;

        pdfName = fullfile(output_folder, ...
            sprintf('%s_%s_arrowup_arrowdown_ttest_CBPT_spectrogram.pdf', ...
            ptID, safeAreaName));

        print(f, pdfName, '-dpdf', '-r300', '-painters');

        close(f);

    end

end

function dataOut = smooth_time_freq_matrix(dataIn, timeSigma, freqSigma, doTime, doFreq)

    dataOut = dataIn;

    if doTime && timeSigma > 0
        dataOut = smooth_time_only_matrix(dataOut, timeSigma);
    end

    if doFreq && freqSigma > 0
        dataOut = smooth_freq_only_matrix(dataOut, freqSigma);
    end

end

function dataOut = smooth_time_only_matrix(dataIn, sigma)

    if sigma <= 0
        dataOut = dataIn;
        return;
    end

    halfWidth = ceil(4 * sigma);
    x = -halfWidth:halfWidth;

    kernel = exp(-(x.^2) ./ (2 * sigma^2));
    kernel = kernel ./ sum(kernel);

    dataOut = nan(size(dataIn));

    for fi = 1:size(dataIn, 1)

        y = dataIn(fi, :);

        nanMask = isnan(y);

        yFilled = y;
        yFilled(nanMask) = 0;

        weight = double(~nanMask);

        smoothY = conv(yFilled, kernel, 'same');
        smoothW = conv(weight, kernel, 'same');

        valid = smoothW > 0;

        smoothY(valid) = smoothY(valid) ./ smoothW(valid);
        smoothY(~valid) = NaN;

        dataOut(fi, :) = smoothY;

    end

end

function dataOut = smooth_freq_only_matrix(dataIn, sigma)

    if sigma <= 0
        dataOut = dataIn;
        return;
    end

    halfWidth = ceil(4 * sigma);
    x = -halfWidth:halfWidth;

    kernel = exp(-(x.^2) ./ (2 * sigma^2));
    kernel = kernel ./ sum(kernel);

    dataOut = nan(size(dataIn));

    for ti = 1:size(dataIn, 2)

        y = dataIn(:, ti);

        nanMask = isnan(y);

        yFilled = y;
        yFilled(nanMask) = 0;

        weight = double(~nanMask);

        smoothY = conv(yFilled, kernel(:), 'same');
        smoothW = conv(weight, kernel(:), 'same');

        valid = smoothW > 0;

        smoothY(valid) = smoothY(valid) ./ smoothW(valid);
        smoothY(~valid) = NaN;

        dataOut(:, ti) = smoothY;

    end

end

function [sigMask, ...
          sigMaskPos, ...
          sigMaskNeg, ...
          pMap, ...
          tMap, ...
          clusterThresholdPos, ...
          clusterThresholdNeg, ...
          realClusterTable, ...
          maxNullPos, ...
          maxNullNeg] = run_ttest_cbpt_sign_separated( ...
            areaTrialData, ...
            cond1Idx, ...
            cond2Idx, ...
            freqVecValid, ...
            timeVecValid, ...
            nPerm, ...
            alphaBin, ...
            alphaCluster, ...
            minClusterSize, ...
            minClusterDurationMs)

    nFreq = size(areaTrialData, 1);
    nTime = size(areaTrialData, 2);
    nTrials = size(areaTrialData, 3);

    fprintf('Fitting real t-tests...\n');

    [tMap, pMap] = compute_ttest_maps_fast( ...
        areaTrialData, ...
        cond1Idx, ...
        cond2Idx);

    candidatePos = pMap < alphaBin & tMap > 0;
    candidateNeg = pMap < alphaBin & tMap < 0;

    candidatePos(isnan(candidatePos)) = false;
    candidateNeg(isnan(candidateNeg)) = false;

    fprintf('Positive candidate bins p < %.3f: %d\n', ...
        alphaBin, sum(candidatePos(:)));

    fprintf('Negative candidate bins p < %.3f: %d\n', ...
        alphaBin, sum(candidateNeg(:)));

    fprintf('Minimum p-value in map: %.6f\n', ...
        min(pMap(:), [], 'omitnan'));

    [realStatsPos, realPixPos] = get_cluster_stats(candidatePos, tMap, minClusterSize, "positive");
    [realStatsNeg, realPixNeg] = get_cluster_stats(candidateNeg, tMap, minClusterSize, "negative");

    fprintf('Positive candidate clusters after min size %d: %d\n', ...
        minClusterSize, numel(realStatsPos));

    fprintf('Negative candidate clusters after min size %d: %d\n', ...
        minClusterSize, numel(realStatsNeg));

    labels = zeros(nTrials,1);

    labels(cond1Idx) = 1;
    labels(cond2Idx) = 2;

    validTrials = find(labels > 0);
    labelsValid = labels(validTrials);

    maxNullPos = zeros(nPerm,1);
    maxNullNeg = zeros(nPerm,1);

    fprintf('Running %d permutations...\n', nPerm);

    for permIdx = 1:nPerm

        shuffledLabels = labelsValid(randperm(numel(labelsValid)));

        permCond1Idx = validTrials(shuffledLabels == 1);
        permCond2Idx = validTrials(shuffledLabels == 2);

        [permTMap, permPMap] = compute_ttest_maps_fast( ...
            areaTrialData, ...
            permCond1Idx, ...
            permCond2Idx);

        permCandidatePos = permPMap < alphaBin & permTMap > 0;
        permCandidateNeg = permPMap < alphaBin & permTMap < 0;

        permCandidatePos(isnan(permCandidatePos)) = false;
        permCandidateNeg(isnan(permCandidateNeg)) = false;

        [permStatsPos, ~] = get_cluster_stats(permCandidatePos, permTMap, minClusterSize, "positive");
        [permStatsNeg, ~] = get_cluster_stats(permCandidateNeg, permTMap, minClusterSize, "negative");

        if isempty(permStatsPos)
            maxNullPos(permIdx) = 0;
        else
            maxNullPos(permIdx) = max(permStatsPos);
        end

        if isempty(permStatsNeg)
            maxNullNeg(permIdx) = 0;
        else
            maxNullNeg(permIdx) = max(permStatsNeg);
        end

        if mod(permIdx,100) == 0
            fprintf('Permutation %d / %d complete\n', ...
                permIdx, nPerm);
        end

    end

    clusterThresholdPos = prctile(maxNullPos, 100*(1-alphaCluster));
    clusterThresholdNeg = prctile(maxNullNeg, 100*(1-alphaCluster));

    sigMaskPos = false(nFreq, nTime);
    sigMaskNeg = false(nFreq, nTime);

    clusterRows = {};

    nSigPos = 0;
    nSigNeg = 0;

    for ci = 1:numel(realStatsPos)

        pix = realPixPos{ci};

        [freqRange, timeRange] = get_cluster_ranges(pix, nFreq, freqVecValid, timeVecValid);

        clusterDurationMs = timeRange(2) - timeRange(1) + 1;

        isSig = realStatsPos(ci) > clusterThresholdPos && clusterDurationMs >= minClusterDurationMs;

        if isSig

            nSigPos = nSigPos + 1;

            sigMaskPos(pix) = true;

            clusterRows(end+1, :) = { ...
                "positive_arrowup_greater_than_arrowdown", ...
                ci, ...
                realStatsPos(ci), ...
                clusterThresholdPos, ...
                numel(pix), ...
                clusterDurationMs, ...
                freqRange(1), ...
                freqRange(2), ...
                timeRange(1), ...
                timeRange(2), ...
                isSig};

        end

    end

    for ci = 1:numel(realStatsNeg)

        pix = realPixNeg{ci};

        [freqRange, timeRange] = get_cluster_ranges(pix, nFreq, freqVecValid, timeVecValid);

        clusterDurationMs = timeRange(2) - timeRange(1) + 1;

        isSig = realStatsNeg(ci) > clusterThresholdNeg && clusterDurationMs >= minClusterDurationMs;

        if isSig

            nSigNeg = nSigNeg + 1;

            sigMaskNeg(pix) = true;

            clusterRows(end+1, :) = { ...
                "negative_arrowdown_greater_than_arrowup", ...
                ci, ...
                realStatsNeg(ci), ...
                clusterThresholdNeg, ...
                numel(pix), ...
                clusterDurationMs, ...
                freqRange(1), ...
                freqRange(2), ...
                timeRange(1), ...
                timeRange(2), ...
                isSig};

        end

    end

    sigMask = sigMaskPos | sigMaskNeg;

    if isempty(clusterRows)

        realClusterTable = table();

    else

        realClusterTable = cell2table(clusterRows, ...
            'VariableNames', { ...
            'direction', ...
            'clusterIndex', ...
            'clusterMass', ...
            'clusterThreshold', ...
            'clusterSizeBins', ...
            'clusterDurationMs', ...
            'freqMin', ...
            'freqMax', ...
            'timeMin', ...
            'timeMax', ...
            'isSignificant'});

    end

    fprintf('Positive cluster threshold = %.4f\n', clusterThresholdPos);
    fprintf('Negative cluster threshold = %.4f\n', clusterThresholdNeg);

    fprintf('Significant positive clusters = %d / %d\n', ...
        nSigPos, numel(realStatsPos));

    fprintf('Significant negative clusters = %d / %d\n', ...
        nSigNeg, numel(realStatsNeg));

end

function [clusterStats, clusterPixels] = get_cluster_stats(candidateMask, tMap, minClusterSize, direction)

    CC = bwconncomp(candidateMask, 8);

    clusterStats = [];
    clusterPixels = {};

    for ci = 1:CC.NumObjects

        pix = CC.PixelIdxList{ci};

        if numel(pix) < minClusterSize
            continue;
        end

        switch direction

            case "positive"
                vals = tMap(pix);
                clusterStat = sum(vals(vals > 0), 'omitnan');

            case "negative"
                vals = tMap(pix);
                clusterStat = sum(abs(vals(vals < 0)), 'omitnan');

            otherwise
                vals = tMap(pix);
                clusterStat = sum(abs(vals), 'omitnan');

        end

        if ~isnan(clusterStat) && clusterStat > 0

            clusterStats(end+1,1) = clusterStat;
            clusterPixels{end+1,1} = pix;

        end

    end

end

function [freqRange, timeRange] = get_cluster_ranges(pix, nFreq, freqVecValid, timeVecValid)

    [freqSub, timeSub] = ind2sub([nFreq, numel(timeVecValid)], pix);

    freqRange = [min(freqVecValid(freqSub)), max(freqVecValid(freqSub))];
    timeRange = [min(timeVecValid(timeSub)), max(timeVecValid(timeSub))];

end

function [tMap, pMap] = compute_ttest_maps_fast( ...
    areaTrialData, ...
    cond1Idx, ...
    cond2Idx)

    cond1Data = areaTrialData(:,:,cond1Idx);
    cond2Data = areaTrialData(:,:,cond2Idx);

    nCond1 = sum(~isnan(cond1Data), 3);
    nCond2 = sum(~isnan(cond2Data), 3);

    meanCond1 = mean(cond1Data, 3, 'omitnan');
    meanCond2 = mean(cond2Data, 3, 'omitnan');

    varCond1 = var(cond1Data, 0, 3, 'omitnan');
    varCond2 = var(cond2Data, 0, 3, 'omitnan');

    se = sqrt((varCond1 ./ nCond1) + (varCond2 ./ nCond2));

    tMap = (meanCond1 - meanCond2) ./ se;

    dfNumer = ((varCond1 ./ nCond1) + (varCond2 ./ nCond2)).^2;
    dfDenom = ((varCond1 ./ nCond1).^2 ./ (nCond1 - 1)) + ((varCond2 ./ nCond2).^2 ./ (nCond2 - 1));

    df = dfNumer ./ dfDenom;

    pMap = 2 .* tcdf(-abs(tMap), df);

    bad = nCond1 < 2 | nCond2 < 2 | se == 0 | isnan(se) | isnan(df) | df <= 0;

    tMap(bad) = NaN;
    pMap(bad) = NaN;

end

function cleanMask = remove_short_time_runs(sigMask, timeVecValid, minRunMs)

    cleanMask = false(size(sigMask));

    for fi = 1:size(sigMask, 1)

        rowMask = sigMask(fi, :);

        CC = bwconncomp(rowMask, 4);

        for ci = 1:CC.NumObjects

            pix = CC.PixelIdxList{ci};

            runTime = timeVecValid(pix);

            runDurationMs = max(runTime) - min(runTime) + 1;

            if runDurationMs >= minRunMs
                cleanMask(fi, pix) = true;
            end

        end

    end

end
