% choice based spectrogram visualization with t-test CBPT
% statistics: average channels within brain area per trial, then arrowup vs arrowdown t-test
% visualization: arrowup/arrowdown average spectrograms, nonsignificant TF bins dimmed

clc;
clear;
close all;

output_folder = fullfile('\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_12_CBPT_choice\');
input_folder  = fullfile('\\155.100.91.44\d\Data\Nill\starling\spectrograms\choice\');

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

    choice = string(bhvData.choice);

    arrowUpTrials   = find(choice == "arrowup");
    arrowDownTrials = find(choice == "arrowdown");

    condNames  = {'arrowup', 'arrowdown'};
    condTrials = {arrowUpTrials, arrowDownTrials};

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

    areaCondData = cell(nAreas, 2);
    areaSigMask = cell(nAreas, 1);
    areaPMap = cell(nAreas, 1);
    areaTMap = cell(nAreas, 1);
    rowClimVals = cell(nAreas, 1);

    for a = 1:nAreas

        areaName = brainAreas(a);
        areaChans = find(cleanAnat == areaName);

        fprintf('\nRunning t-test CBPT for area: %s | channels: %d\n', areaName, numel(areaChans));

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

                chanData(:, :, k) = tmp;

            end

            areaTrialData(:, :, tr) = mean(chanData, 3, 'omitnan');

        end

        arrowUpLocalIdx   = find(statsChoice == "arrowup");
        arrowDownLocalIdx = find(statsChoice == "arrowdown");

        if isempty(arrowUpLocalIdx) || isempty(arrowDownLocalIdx)

            fprintf('Skipping %s because one condition is empty.\n', areaName);

            sigMaskSmall = false(numel(validFreqIdx), numel(validTimeIdx));
            pMapSmall = nan(numel(validFreqIdx), numel(validTimeIdx));
            tMapSmall = nan(numel(validFreqIdx), numel(validTimeIdx));

        else

            [sigMaskSmall, pMapSmall, tMapSmall] = run_ttest_cbpt( ...
                areaTrialData, arrowUpLocalIdx, arrowDownLocalIdx, ...
                nPerm, alphaBin, alphaCluster, minClusterSize);

        end

        fullSigMask = false(nFreq, nTime);
        fullPMap = nan(nFreq, nTime);
        fullTMap = nan(nFreq, nTime);

        fullSigMask(validFreqIdx, validTimeIdx) = sigMaskSmall;
        fullPMap(validFreqIdx, validTimeIdx) = pMapSmall;
        fullTMap(validFreqIdx, validTimeIdx) = tMapSmall;

        areaSigMask{a} = fullSigMask;
        areaPMap{a} = fullPMap;
        areaTMap{a} = fullTMap;

        rowVals = [];

        for c = 1:2

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

        save(fullfile(output_folder, sprintf('%s_%s_choice_ttest_CBPT_results.mat', ptID, safeAreaName)), ...
            'sigMaskSmall', 'fullSigMask', ...
            'pMapSmall', 'tMapSmall', 'fullPMap', 'fullTMap', ...
            'freqVec', 'validFreqIdx', 'validTimeIdx', ...
            'areaName', 'areaChans', ...
            'nPerm', 'alphaBin', 'alphaCluster', 'minClusterSize', '-v7.3');

    end

    rowsPerPage = 5;
    nPages = ceil(nAreas / rowsPerPage);

    for pageNum = 1:nPages

        areaStart = (pageNum - 1) * rowsPerPage + 1;
        areaEnd = min(pageNum * rowsPerPage, nAreas);
        pageAreas = areaStart:areaEnd;
        nPageRows = numel(pageAreas);

        f = figure('Visible', 'off', ...
            'Position', [100 100 1400 400 * nPageRows]);

        t = tiledlayout(nPageRows, 2, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        title(t, sprintf('%s | Choice Spectrogram with t-test CBPT | Page %d', ptID, pageNum), ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');

        for rr = 1:nPageRows

            a = pageAreas(rr);

            for c = 1:2

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
                xline(1000, 'r', 'LineWidth', 0.5);
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

                if c == 2
                    h = colorbar(ax, 'eastoutside');
                    h.Label.String = 'power';
                    h.FontSize = 8;
                end

            end

        end

        set(f, 'Renderer', 'opengl');

        exportgraphics(f, ...
            fullfile(output_folder, sprintf('%s_arrowup_arrowdown_ttest_CBPT_spectrogram_page_%02d.pdf', ptID, pageNum)), ...
            'ContentType', 'image', ...
            'BackgroundColor', 'white', ...
            'Resolution', 600);

        close(f);

    end

end

function [sigMask, pMap, tMap] = run_ttest_cbpt(areaTrialData, cond1Idx, cond2Idx, nPerm, alphaBin, alphaCluster, minClusterSize)

    nFreq = size(areaTrialData, 1);
    nTime = size(areaTrialData, 2);
    nTrials = size(areaTrialData, 3);

    fprintf('Fitting real t-tests...\n');

    [tMap, pMap] = compute_ttest_maps(areaTrialData, cond1Idx, cond2Idx);

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

        clusterStat = sum(abs(tMap(pix)), 'omitnan');

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

    labels = zeros(nTrials, 1);
    labels(cond1Idx) = 1;
    labels(cond2Idx) = 2;

    validTrials = find(labels > 0);
    labelsValid = labels(validTrials);

    maxNull = zeros(nPerm, 1);

    fprintf('Running %d permutations...\n', nPerm);

    for permIdx = 1:nPerm

        shuffledLabels = labelsValid(randperm(numel(labelsValid)));

        permCond1Idx = validTrials(shuffledLabels == 1);
        permCond2Idx = validTrials(shuffledLabels == 2);

        [permTMap, permPMap] = compute_ttest_maps(areaTrialData, permCond1Idx, permCond2Idx);

        permCandidateMask = permPMap < alphaBin;
        permCandidateMask(isnan(permCandidateMask)) = false;

        permCC = bwconncomp(permCandidateMask, 8);

        permClusterStats = [];

        for ci = 1:permCC.NumObjects

            pix = permCC.PixelIdxList{ci};

            if numel(pix) < minClusterSize
                continue;
            end

            clusterStat = sum(abs(permTMap(pix)), 'omitnan');

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

function [tMap, pMap] = compute_ttest_maps(areaTrialData, cond1Idx, cond2Idx)

    nFreq = size(areaTrialData, 1);
    nTime = size(areaTrialData, 2);

    tMap = nan(nFreq, nTime);
    pMap = nan(nFreq, nTime);

    cond1Data = areaTrialData(:, :, cond1Idx);
    cond2Data = areaTrialData(:, :, cond2Idx);

    for fi = 1:nFreq

        for ti = 1:nTime

            x = squeeze(cond1Data(fi, ti, :));
            y = squeeze(cond2Data(fi, ti, :));

            x = x(~isnan(x));
            y = y(~isnan(y));

            if numel(x) < 2 || numel(y) < 2
                continue;
            end

            try
                [~, p, ~, stats] = ttest2(x, y, 'Vartype', 'unequal');

                pMap(fi, ti) = p;
                tMap(fi, ti) = stats.tstat;

            catch
                pMap(fi, ti) = NaN;
                tMap(fi, ti) = NaN;
            end

        end

    end

end