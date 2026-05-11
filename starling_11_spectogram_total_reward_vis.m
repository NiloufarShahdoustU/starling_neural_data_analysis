% outcome based spectrogram visualization
% One PDF per participant/page
% rows = brain areas, with left/right combined
% columns = win / lose
% each row has its own colorbar shared between win and lose
% saves everything directly in one output folder

clc;
clear;
close all;

output_folder = fullfile('\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_11_spectogram_total_reward_vis\');
input_folder  = fullfile('\\155.100.91.44\d\Data\Nill\starling\spectrograms\total_reward\');

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

files = dir(fullfile(input_folder, '*_spectrogram_data.mat'));

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
    fWin = [1 200];

    outcome = string(bhvData.outcome);

    winTrials  = find(outcome == "win");
    loseTrials = find(outcome == "lose");

    condNames  = {'win', 'lose'};
    condTrials = {winTrials, loseTrials};

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
    rowClimVals = cell(nAreas, 1);

    for a = 1:nAreas

        areaName = brainAreas(a);
        areaChans = find(cleanAnat == areaName);

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

        title(t, sprintf('%s | Outcome Spectrogram by Brain Area | Page %d', ptID, pageNum), ...
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

                % remove frequencies below 5 Hz
                freqMask = freqVec >= 5;
                
                % remove last 100 ms from time
                timeVec = 1:size(dataToPlot,2);
                timeVec = timeVec(1:end-100);
                
                dataPlot = dataToPlot(freqMask, 1:end-100);
                freqPlot = freqVec(freqMask);
                
                imagesc(timeVec, freqPlot, dataPlot);
                
                set(gca, 'YDir', 'normal', 'YScale', 'log');
                
                ylim([5, fWin(2)]);

                if ~isempty(rowClimVals{a})
                    caxis(rowClimVals{a});
                end

                axis tight;

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
                    xlabel('Time');
                end

                if c == 2
                    h = colorbar(ax, 'eastoutside');
                    h.Label.String = 'power';
                    h.FontSize = 8;
                end

            end

        end

        set(f, 'Renderer', 'painters');

        exportgraphics(f, ...
            fullfile(output_folder, sprintf('%s_win_lose_spectrogram_page_%02d.pdf', ptID, pageNum)), ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'none', ...
            'Resolution', 600);

        close(f);

    end

end