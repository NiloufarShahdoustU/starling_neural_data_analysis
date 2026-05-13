% flip based spectrogram visualization with one-way ANOVA CBPT
% conditions: uniform / low / high
% statistics: average channels within brain area per trial, then one-way ANOVA per TF bin
% visualization: condition averages, nonsignificant TF bins dimmed

clc;
clear;
close all;

output_folder = fullfile('\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_12_CBPT_flip_v2\');
input_folder  = fullfile('\\155.100.91.44\d\Data\Nill\starling\spectrograms\flip\');

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

    cue = lower(string(bhvData.distribution));

    uniformTrials = find(cue == "uniform");
    lowTrials     = find(cue == "low");
    highTrials    = find(cue == "high");

    condNames  = {'uniform', 'low', 'high'};
    condTrials = {uniformTrials, lowTrials, highTrials};

    cleanAnat = strings(nChans, 1);

    for ch = 1:nChans

        thisLoc = string(anatomicalLocs{ch});

        thisLoc = regexprep(thisLoc, 'left', '', 'ignorecase');
        thisLoc = regexprep(thisLoc, 'right', '', 'ignorecase');
        thisLoc = regexprep(thisLoc, '\blh\b', '', 'ignorecase');
        thisLoc = regexprep(thisLoc, '\brh\b', '', 'ignorecase');
        thisLoc = regexprep(thisLoc, '(^|\W)L(\W|$)', ' ', 'ignorecase');
        thisLoc = regexprep(thisLoc, '(^|\W)R(\W|$)', ' ', 'ignorecase');

        thisLoc = regexprep(thisLoc, '[-_]', ' ');
        thisLoc = strtrim(regexprep(thisLoc, '\s+', ' '));

        cleanAnat(ch) = thisLoc;

    end

    brainAreas = unique(cleanAnat, 'stable');
    brainAreas(brainAreas == "") = [];

    for a = 1:numel(brainAreas)

        areaName = brainAreas(a);
        areaChans = find(cleanAnat == areaName);

        safeAreaName = matlab.lang.makeValidName(char(areaName));

        fprintf('\nRunning one-way ANOVA CBPT for area: %s | channels: %d\n', ...
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

        statsTrials = [uniformTrials; lowTrials; highTrials];

        groupLabels = [
            ones(numel(uniformTrials), 1);
            2 * ones(numel(lowTrials), 1);
            3 * ones(numel(highTrials), 1)
        ];

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

        if numel(unique(groupLabels)) < 3 || ...
                sum(groupLabels == 1) < 2 || ...
                sum(groupLabels == 2) < 2 || ...
                sum(groupLabels == 3) < 2

            fprintf('Skipping %s because one condition has fewer than 2 trials.\n', areaName);

            sigMaskSmall = false(numel(validFreqIdx), numel(validTimeIdx));
            pMapSmall = nan(numel(validFreqIdx), numel(validTimeIdx));
            fMapSmall = nan(numel(validFreqIdx), numel(validTimeIdx));

            clusterThreshold = NaN;
            realClusterTable = table();
            maxNull = [];

        else

            [sigMaskSmall, ...
             pMapSmall, ...
             fMapSmall, ...
             clusterThreshold, ...
             realClusterTable, ...
             maxNull] = run_anova_cbpt( ...
                areaTrialData, ...
                groupLabels, ...
                freqVecValid, ...
                validTimeIdx, ...
                nPerm, ...
                alphaBin, ...
                alphaCluster, ...
                minClusterSize, ...
                minClusterDurationMs);

            sigMaskSmall = remove_short_time_runs(sigMaskSmall, validTimeIdx, minSigPixelRunMs);

        end

        fullSigMask = false(nFreq, nTime);
        fullPMap = nan(nFreq, nTime);
        fullFMap = nan(nFreq, nTime);

        fullSigMask(validFreqIdx, validTimeIdx) = sigMaskSmall;
        fullPMap(validFreqIdx, validTimeIdx) = pMapSmall;
        fullFMap(validFreqIdx, validTimeIdx) = fMapSmall;

        matName = fullfile(output_folder, ...
            sprintf('%s_%s_flip_anova_CBPT_results.mat', ...
            ptID, safeAreaName));

        save(matName, ...
            'sigMaskSmall', ...
            'fullSigMask', ...
            'pMapSmall', ...
            'fMapSmall', ...
            'fullPMap', ...
            'fullFMap', ...
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
            'clusterThreshold', ...
            'realClusterTable', ...
            'maxNull', ...
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

        areaCondDataSmoothed = cell(1, 3);
        rowVals = [];

        for c = 1:3

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
            'Position',[1 1 24 9], ...
            'PaperUnits','inches', ...
            'PaperSize',[24 9], ...
            'PaperPosition',[0 0 24 9]);

        sgtitle(sprintf('%s | %s | time-frequency flip spectrogram ANOVA CBPT', ...
            ptID, areaName), ...
            'FontWeight','bold', ...
            'Interpreter','none');

        axPos = [
            0.05 0.16 0.25 0.70
            0.38 0.16 0.25 0.70
            0.71 0.16 0.25 0.70
        ];

        for c = 1:3

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
            sprintf('%s_%s_uniform_low_high_flip_anova_CBPT_spectrogram.pdf', ...
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
          pMap, ...
          fMap, ...
          clusterThreshold, ...
          realClusterTable, ...
          maxNull] = run_anova_cbpt( ...
            areaTrialData, ...
            groupLabels, ...
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

    fprintf('Fitting real one-way ANOVA maps...\n');

    [fMap, pMap] = compute_anova_maps(areaTrialData, groupLabels);

    candidateMask = pMap < alphaBin;
    candidateMask(isnan(candidateMask)) = false;

    fprintf('Number of candidate bins p < %.3f: %d\n', alphaBin, sum(candidateMask(:)));
    fprintf('Minimum p-value in map: %.6f\n', min(pMap(:), [], 'omitnan'));

    [realClusterStats, realClusterPixels] = get_anova_cluster_stats(candidateMask, fMap, minClusterSize);

    fprintf('Candidate clusters after min size %d: %d\n', minClusterSize, numel(realClusterStats));

    maxNull = zeros(nPerm, 1);

    fprintf('Running %d permutations...\n', nPerm);

    for permIdx = 1:nPerm

        shuffledLabels = groupLabels(randperm(nTrials));

        [permFMap, permPMap] = compute_anova_maps(areaTrialData, shuffledLabels);

        permCandidateMask = permPMap < alphaBin;
        permCandidateMask(isnan(permCandidateMask)) = false;

        [permClusterStats, ~] = get_anova_cluster_stats(permCandidateMask, permFMap, minClusterSize);

        if isempty(permClusterStats)
            maxNull(permIdx) = 0;
        else
            maxNull(permIdx) = max(permClusterStats);
        end

        if mod(permIdx, 100) == 0
            fprintf('Permutation %d / %d complete\n', permIdx, nPerm);
        end

    end

    clusterThreshold = prctile(maxNull, 100 * (1 - alphaCluster));

    sigMask = false(nFreq, nTime);
    clusterRows = {};
    nSig = 0;

    for ci = 1:numel(realClusterStats)

        pix = realClusterPixels{ci};

        [freqRange, timeRange] = get_cluster_ranges(pix, nFreq, freqVecValid, timeVecValid);

        clusterDurationMs = timeRange(2) - timeRange(1) + 1;

        isSig = realClusterStats(ci) > clusterThreshold && clusterDurationMs >= minClusterDurationMs;

        if isSig

            nSig = nSig + 1;

            sigMask(pix) = true;

            clusterRows(end+1, :) = { ...
                "anova_uniform_low_high", ...
                ci, ...
                realClusterStats(ci), ...
                clusterThreshold, ...
                numel(pix), ...
                clusterDurationMs, ...
                freqRange(1), ...
                freqRange(2), ...
                timeRange(1), ...
                timeRange(2), ...
                isSig};

        end

    end

    if isempty(clusterRows)

        realClusterTable = table();

    else

        realClusterTable = cell2table(clusterRows, ...
            'VariableNames', { ...
            'test', ...
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

    fprintf('Cluster threshold = %.4f\n', clusterThreshold);
    fprintf('Significant clusters = %d / %d\n', nSig, numel(realClusterStats));

end

function [clusterStats, clusterPixels] = get_anova_cluster_stats(candidateMask, fMap, minClusterSize)

    CC = bwconncomp(candidateMask, 8);

    clusterStats = [];
    clusterPixels = {};

    for ci = 1:CC.NumObjects

        pix = CC.PixelIdxList{ci};

        if numel(pix) < minClusterSize
            continue;
        end

        clusterStat = sum(fMap(pix), 'omitnan');

        if ~isnan(clusterStat) && clusterStat > 0

            clusterStats(end+1, 1) = clusterStat;
            clusterPixels{end+1, 1} = pix;

        end

    end

end

function [freqRange, timeRange] = get_cluster_ranges(pix, nFreq, freqVecValid, timeVecValid)

    [freqSub, timeSub] = ind2sub([nFreq, numel(timeVecValid)], pix);

    freqRange = [min(freqVecValid(freqSub)), max(freqVecValid(freqSub))];
    timeRange = [min(timeVecValid(timeSub)), max(timeVecValid(timeSub))];

end

function [fMap, pMap] = compute_anova_maps(areaTrialData, groupLabels)

    nFreq = size(areaTrialData, 1);
    nTime = size(areaTrialData, 2);

    fMap = nan(nFreq, nTime);
    pMap = nan(nFreq, nTime);

    groupLabels = groupLabels(:);

    for fi = 1:nFreq

        for ti = 1:nTime

            y = squeeze(areaTrialData(fi, ti, :));
            y = y(:);

            validIdx = ~isnan(y) & ~isnan(groupLabels);

            yValid = y(validIdx);
            gValid = groupLabels(validIdx);

            if numel(unique(gValid)) < 3
                continue;
            end

            if sum(gValid == 1) < 2 || sum(gValid == 2) < 2 || sum(gValid == 3) < 2
                continue;
            end

            groupMeans = [
                mean(yValid(gValid == 1), 'omitnan');
                mean(yValid(gValid == 2), 'omitnan');
                mean(yValid(gValid == 3), 'omitnan')
            ];

            grandMean = mean(yValid, 'omitnan');

            n1 = sum(gValid == 1);
            n2 = sum(gValid == 2);
            n3 = sum(gValid == 3);

            ssBetween = n1 * (groupMeans(1) - grandMean)^2 + ...
                        n2 * (groupMeans(2) - grandMean)^2 + ...
                        n3 * (groupMeans(3) - grandMean)^2;

            ssWithin = sum((yValid(gValid == 1) - groupMeans(1)).^2) + ...
                       sum((yValid(gValid == 2) - groupMeans(2)).^2) + ...
                       sum((yValid(gValid == 3) - groupMeans(3)).^2);

            dfBetween = 3 - 1;
            dfWithin = numel(yValid) - 3;

            if dfWithin > 0 && ssWithin > 0
                fMap(fi, ti) = (ssBetween / dfBetween) / (ssWithin / dfWithin);
                pMap(fi, ti) = 1 - fcdf(fMap(fi, ti), dfBetween, dfWithin);
            end

        end

    end

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
