clc;
clear;
close all;
%%
input_folder = '\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_9_spectogram_4_areas_prep';
output_folder = '\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_9_spectogram_4_areas_visualization';

targetAreas = {'Hippocampus', 'MTG', 'Amygdala', 'SFG'};
fWin = [1 200]; % frequency range of interest

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

areaData = cell(1, numel(targetAreas));
freqRef = []; % to store frequencies once

files = dir(fullfile(input_folder, '*.mat'));

for f = 1:numel(files)
    filePath = fullfile(input_folder, files(f).name);
    fprintf('Loading %s...\n', files(f).name);

    S = load(filePath);

    dataToSave = S.dataToSave;                        % [nCh x nFreq x nTime]
    SelectedAnatomicalLocs = S.SelectedAnatomicalLoc; % {1 x nCh}
    freqToSave = S.freqToSave;                        % [1 x nFreq]

    if isempty(freqRef)
        freqRef = freqToSave; % store reference frequencies
    end

    for a = 1:numel(targetAreas)
        areaName = targetAreas{a};
        idx = find(contains(SelectedAnatomicalLocs, areaName, 'IgnoreCase', true));

        if ~isempty(idx)
            for i = 1:numel(idx)
                thisData = squeeze(dataToSave(idx(i), :, :)); % [nFreq x nTime]
                areaData{a}{end+1} = thisData;
            end
        end
    end
end

%% Compute mean spectrograms
meanSpectra = cell(1, numel(targetAreas));
for a = 1:numel(targetAreas)
    if isempty(areaData{a})
        warning('No data found for %s', targetAreas{a});
        continue;
    end
    stacked = cat(3, areaData{a}{:});        % [nFreq x nTime x nCh_total]
    meanSpectra{a} = mean(stacked, 3, 'omitnan');
end

%% Visualization (individual CLim for each)

figure('Color','w','Position',[200 200 1200 300]);

for a = 1:numel(targetAreas)
    subplot(1, numel(targetAreas), a);

    if isempty(meanSpectra{a})
        title([targetAreas{a} ' (no data)']);
        continue;
    end

    % Compute color limits separately for each area
    vals = meanSpectra{a}(:);
    vals = vals(~isnan(vals));
    lowCut  = prctile(vals, 20);
    highCut = prctile(vals, 87);
    clims = [lowCut, highCut];

    % Plot
    imagesc(1:size(meanSpectra{a}, 2), freqRef, meanSpectra{a});
    set(gca, 'YDir','normal', 'CLim', clims, 'YScale', 'log');
    ylim([fWin(1) fWin(2)]);
    xline(1000, 'r', 'LineWidth', 2);

    % Title and labels
    title(targetAreas{a}, 'Interpreter','none');
    xlabel('time (ms)');

    % Only first subplot has Y-axis label
    if a == 1
        ylabel('frequency (Hz)');
    else
        ylabel('');
    end

    % Log scale ticks: 10^0, 10^1, 10^2
    yticks([1 10 100]);
    yticklabels({'10^0', '10^1', '10^2'});

    % X-axis ticks: 1→−1000, 1000→0, 2000→1000
    xticks([1 1000 2000]);
    xticklabels({'-1000','0','1000'});

    % Colorbar
    cb = colorbar;
    if a == numel(targetAreas)
        cb.Label.String = 'power (a.u.)';
        cb.Label.FontSize = 8;
    else
        cb.Label.String = '';
    end

    axis square;

    fprintf('%s color limits: [%.4f, %.4f]\n', targetAreas{a}, clims(1), clims(2));
end

% ------------------------------------------------------
% Bring subplots closer together
% ------------------------------------------------------
subplotHandles = findall(gcf, 'Type', 'axes');
for i = 1:numel(subplotHandles)
    pos = get(subplotHandles(i), 'Position');
    pos(3) = pos(3) * 1.15;   % increase width
    pos(1) = pos(1) - 0.025;  % shift slightly left
    set(subplotHandles(i), 'Position', pos);
end

sgtitle('mean spectrograms');
set(gcf, 'Renderer', 'painters');

pdf_path = fullfile(output_folder, 'mean_spectrograms.pdf');
exportgraphics(gcf, pdf_path, 'ContentType', 'vector', 'Resolution', 300);
