% cue based spectrogram visualization with one-way ANOVA CBPT
% conditions: uniform / low / high
% statistics: average channels within brain area per trial, then one-way ANOVA per TF bin
% visualization: condition averages, nonsignificant TF bins dimmed

clc;
clear;
close all;

output_folder = fullfile('\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_12_CBPT_cue\');
input_folder  = fullfile('\\155.100.91.44\d\Data\Nill\starling\spectrograms\cue\');

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

files = dir(fullfile(input_folder, '*_spectrogram_data.mat'));

nPerm = 1000;

alphaBin = 0.05;
alphaCluster = 0.05;
minClusterSize = 5;

fWin = [5 200];
removeLastMs = 100;

rng(1);

for p = 1:numel(files)

    file_name = files(p).name;
    file_path = fullfile(files(p).folder, file_name);

    ptID = erase(file_name, '_spectrogram_data.mat');

    fprintf('\n--- Processing ptID: %s ---\n', ptID);

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

    nAreas = numel(brainAreas);

    areaCondData = cell(nAreas, 3);
    areaSigMask = cell(nAreas, 1);
    areaPMap = cell(nAreas, 1);
    areaFMap = cell(nAreas, 1);
    rowClimVals = cell(nAreas, 1);

    for a = 1:nAreas

        areaName = brainAreas(a);
        areaChans = find(cleanAnat == areaName);

        fprintf('\nRunning one-way ANOVA CBPT for area: %s | channels: %d\n', areaName, numel(areaChans));

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

                chanData(:, :, k) = tmp;

            end

            areaTrialData(:, :, tr) = mean(chanData, 3, 'omitnan');

        end

        if numel(unique(groupLabels)) < 3 || ...
                sum(groupLabels == 1) < 2 || ...
                sum(groupLabels == 2) < 2 || ...
                sum(groupLabels == 3) < 2

            fprintf('Skipping %s because one condition has fewer than 2 trials.\n', areaName);

            sigMaskSmall = false(numel(validFreqIdx), numel(validTimeIdx));
            pMapSmall = nan(numel(validFreqIdx), numel(validTimeIdx));
            fMapSmall = nan(numel(validFreqIdx), numel(validTimeIdx));

        else

            [sigMaskSmall, pMapSmall, fMapSmall] = run_anova_cbpt( ...
                areaTrialData, groupLabels, nPerm, alphaBin, alphaCluster, minClusterSize);

        end

        fullSigMask = false(nFreq, nTime);
        fullPMap = nan(nFreq, nTime);
        fullFMap = nan(nFreq, nTime);

        fullSigMask(validFreqIdx, validTimeIdx) = sigMaskSmall;
        fullPMap(validFreqIdx, validTimeIdx) = pMapSmall;
        fullFMap(validFreqIdx, validTimeIdx) = fMapSmall;

        areaSigMask{a} = fullSigMask;
        areaPMap{a} = fullPMap;
        areaFMap{a} = fullFMap;

        rowVals = [];

        for c = 1:3

            trialsIdx = condTrials{c};

            if isempty(trialsIdx)
                areaCondData{a, c} = [];
                continue;
            end

            chanData = [];

            for k = 1:numel(areaChans)

                ch = areaChans(k);

                tmp = mean(Sbl_all{ch}(:, :, trialsIdx), 3, 'omitnan');

                chanData(:, :, k) = tmp;

            end

            areaMean = mean(chanData, 3, 'omitnan');

            areaCondData{a, c} = areaMean;

            rowVals = [rowVals; areaMean(:)];

        end

        rowVals = rowVals(~isnan(rowVals));

        if isempty(rowVals)
            rowClimVals{a} = [];
        else
            rowClimVals{a} = [prctile(rowVals, 5), prctile(rowVals, 95)];
        end

        safeAreaName = matlab.lang.makeValidName(char(areaName));

        save(fullfile(output_folder, sprintf('%s_%s_cue_anova_CBPT_results.mat', ptID, safeAreaName)), ...
            'sigMaskSmall', 'fullSigMask', ...
            'pMapSmall', 'fMapSmall', 'fullPMap', 'fullFMap', ...
            'freqVec', 'validFreqIdx', 'validTimeIdx', ...
            'areaName', 'areaChans', ...
            'condNames', 'nPerm', 'alphaBin', 'alphaCluster', 'minClusterSize', '-v7.3');

    end

    rowsPerPage = 4;
    nPages = ceil(nAreas / rowsPerPage);

    for pageNum = 1:nPages

        areaStart = (pageNum - 1) * rowsPerPage + 1;
        areaEnd = min(pageNum * rowsPerPage, nAreas);
        pageAreas = areaStart:areaEnd;
        nPageRows = numel(pageAreas);

        f = figure('Visible', 'off', ...
            'Position', [100 100 1800 400 * nPageRows]);

        t = tiledlayout(nPageRows, 3, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        title(t, sprintf('%s | Cue Spectrogram with one-way ANOVA CBPT | Page %d', ptID, pageNum), ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');

        for rr = 1:nPageRows

            a = pageAreas(rr);

            for c = 1:3

                ax = nexttile;

                dataToPlot = areaCondData{a, c};

                if isempty(dataToPlot)
                    axis off;
                    title(sprintf('%s: no trials', condNames{c}), ...
                        'Interpreter', 'none', ...
                        'FontSize', 8);
                    continue;
                end

                areaChans = find(cleanAnat == brainAreas(a));
                firstCh = areaChans(1);

                if iscell(freq_all)
                    freqVec = freq_all{firstCh};
                else
                    freqVec = freq_all;
                end

                freqMask = freqVec >= fWin(1) & freqVec <= fWin(2);

                timeVec = 1:size(dataToPlot, 2);
                timeVec = timeVec(1:end-removeLastMs);

                dataPlot = dataToPlot(freqMask, 1:end-removeLastMs);
                freqPlot = freqVec(freqMask);

                sigMaskPlot = areaSigMask{a}(freqMask, 1:end-removeLastMs);

                alphaMask = 0.15 * ones(size(sigMaskPlot));
                alphaMask(sigMaskPlot) = 1.0;

                hImg = imagesc(timeVec, freqPlot, dataPlot);
                set(hImg, 'AlphaData', alphaMask);

                set(gca, 'YDir', 'normal', 'YScale', 'log');

                ylim([fWin(1), fWin(2)]);
                axis tight;

                if ~isempty(rowClimVals{a})
                    caxis(rowClimVals{a});
                end

                hold on;
                xline(500, 'r', 'LineWidth', 0.5);
                hold off;

                set(gca, 'FontSize', 9);

                if c == 1
                    ylabel(char(brainAreas(a)), ...
                        'Interpreter', 'none', ...
                        'FontSize', 9);
                end

                if rr == 1
                    title(condNames{c}, ...
                        'FontWeight', 'bold', ...
                        'FontSize', 12);
                end

                if rr == nPageRows
                    xlabel('Time (ms)');
                end

                if c == 3
                    h = colorbar(ax, 'eastoutside');
                    h.Label.String = 'power';
                    h.FontSize = 8;
                end

            end

        end

        set(f, 'Renderer', 'opengl');

        exportgraphics(f, ...
            fullfile(output_folder, sprintf('%s_uniform_low_high_anova_CBPT_spectrogram_page_%02d.pdf', ptID, pageNum)), ...
            'ContentType', 'image', ...
            'BackgroundColor', 'white', ...
            'Resolution', 600);

        close(f);

    end

end

function [sigMask, pMap, fMap] = run_anova_cbpt(areaTrialData, groupLabels, nPerm, alphaBin, alphaCluster, minClusterSize)

    nFreq = size(areaTrialData, 1);
    nTime = size(areaTrialData, 2);
    nTrials = size(areaTrialData, 3);

    fprintf('Fitting real one-way ANOVA maps...\n');

    [fMap, pMap] = compute_anova_maps(areaTrialData, groupLabels);

    candidateMask = pMap < alphaBin;
    candidateMask(isnan(candidateMask)) = false;

    fprintf('Number of candidate bins p < %.3f: %d\n', alphaBin, sum(candidateMask(:)));
    fprintf('Minimum p-value in map: %.6f\n', min(pMap(:), [], 'omitnan'));

    CC = bwconncomp(candidateMask, 8);

    realClusterStats = [];
    realClusterPixels = {};

    for ci = 1:CC.NumObjects

        pix = CC.PixelIdxList{ci};

        if numel(pix) < minClusterSize
            continue;
        end

        clusterStat = sum(fMap(pix), 'omitnan');

        if ~isnan(clusterStat)
            realClusterStats(end+1, 1) = clusterStat;
            realClusterPixels{end+1, 1} = pix;
        end

    end

    fprintf('Candidate clusters after min size %d: %d\n', minClusterSize, numel(realClusterStats));

    if isempty(realClusterStats)
        sigMask = false(nFreq, nTime);
        fprintf('No candidate clusters found after cluster-size filtering.\n');
        return;
    end

    maxNull = zeros(nPerm, 1);

    fprintf('Running %d permutations...\n', nPerm);

    for permIdx = 1:nPerm

        shuffledLabels = groupLabels(randperm(nTrials));

        [permFMap, permPMap] = compute_anova_maps(areaTrialData, shuffledLabels);

        permCandidateMask = permPMap < alphaBin;
        permCandidateMask(isnan(permCandidateMask)) = false;

        permCC = bwconncomp(permCandidateMask, 8);

        permClusterStats = [];

        for ci = 1:permCC.NumObjects

            pix = permCC.PixelIdxList{ci};

            if numel(pix) < minClusterSize
                continue;
            end

            clusterStat = sum(permFMap(pix), 'omitnan');

            if ~isnan(clusterStat)
                permClusterStats(end+1, 1) = clusterStat;
            end

        end

        if isempty(permClusterStats)
            maxNull(permIdx) = 0;
        else
            maxNull(permIdx) = max(permClusterStats);
        end

        if mod(permIdx, 50) == 0
            fprintf('Permutation %d / %d complete\n', permIdx, nPerm);
        end

    end

    clusterThreshold = prctile(maxNull, 100 * (1 - alphaCluster));

    sigMask = false(nFreq, nTime);

    for ci = 1:numel(realClusterStats)

        if realClusterStats(ci) > clusterThreshold
            sigMask(realClusterPixels{ci}) = true;
        end

    end

    fprintf('Cluster threshold = %.4f\n', clusterThreshold);
    fprintf('Significant clusters = %d / %d\n', sum(realClusterStats > clusterThreshold), numel(realClusterStats));

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

            try
                p = anova1(yValid, gValid, 'off');

                pMap(fi, ti) = p;

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
                end

            catch
                pMap(fi, ti) = NaN;
                fMap(fi, ti) = NaN;
            end

        end

    end

end