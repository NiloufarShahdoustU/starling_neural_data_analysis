clc;
clear;
close all

main_folder = '\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_8_encoding_choice_v2';

d = dir(main_folder);
isub = [d(:).isdir];
subFolders = {d(isub).name}';
subFolders(ismember(subFolders, {'.','..'})) = [];

ptIDs = string(subFolders);
fprintf('Found %d participant folders:\n', numel(ptIDs));
disp(ptIDs);

database = struct();

for i = 1:numel(ptIDs)
    ptID = ptIDs(i);
    fprintf('Now processing %s...\n', ptID);

    mat_path = fullfile(main_folder, ptID, 'encodingResults.mat');
    
    data = load(mat_path);
    fprintf('Loaded %s\n', mat_path);

    fieldname = ['pt' char(ptID)];
    database.(fieldname) = data;
end

%% aggregate cvPseudoR2 by brain area across all participants
models = {'greedy', 'softmax', 'rs'};
brain_data = struct();

ptFields = fieldnames(database);

for i = 1:numel(ptFields)
    ptName = ptFields{i};
    data = database.(ptName);

    if ~isfield(data, 'results') || ~isfield(data, 'NewAnatomicalLocs')
        warning('⚠️ Skipping %s: Missing results/NewAnatomicalLocs', ptName);
        continue;
    end

    locs = string(data.NewAnatomicalLocs);
    results = data.results;

    % Clean region names (replace illegal characters)
    cleanLocs = regexprep(locs, '[^a-zA-Z0-9]', '_');

    % Avoid duplicate regions within a participant
    uniqueLocs = unique(cleanLocs);

    for u = 1:numel(uniqueLocs)
        region = uniqueLocs(u);
        regionName = char(region);              % convert to char
        regionName = matlab.lang.makeValidName(regionName);  % ensure valid field name

        idx = find(cleanLocs == region);

        for m = 1:numel(models)
            mdl = models{m};

            if isfield(results, mdl) && isfield(results.(mdl), 'cvPseudoR2')
                vals = results.(mdl).cvPseudoR2;
                if isempty(vals), continue; end

                validIdx = idx(idx <= numel(vals));
                validVals = vals(validIdx);
                validVals = validVals(~isnan(validVals));

                if ~isempty(validVals)
                    meanVal = mean(validVals);
                    if ~isfield(brain_data, regionName)
                        brain_data.(regionName) = struct();
                    end
                    if ~isfield(brain_data.(regionName), mdl)
                        brain_data.(regionName).(mdl) = [];
                    end
                    brain_data.(regionName).(mdl)(end+1) = meanVal;
                end
            end
        end
    end
end




% --- Combine left/right homologous regions ---
regionNames = fieldnames(brain_data);
combined_brain_data = struct();

for i = 1:numel(regionNames)
    name = regionNames{i};
    % remove left/right prefixes if present
    baseName = regexprep(name, '^(Left_|Right_)', '');
    baseName = strrep(baseName, '__', '_'); % cleanup double underscores

    % ensure valid field name
    baseName = matlab.lang.makeValidName(baseName);

    % initialize if not existing
    if ~isfield(combined_brain_data, baseName)
        combined_brain_data.(baseName) = struct();
    end

    % merge all models' data
    modelFields = fieldnames(brain_data.(name));
    for m = 1:numel(modelFields)
        mdl = modelFields{m};
        if ~isfield(combined_brain_data.(baseName), mdl)
            combined_brain_data.(baseName).(mdl) = [];
        end
        combined_brain_data.(baseName).(mdl) = ...
            [combined_brain_data.(baseName).(mdl), brain_data.(name).(mdl)];
    end
end

% replace original brain_data with combined version
brain_data = combined_brain_data;

%% vis
% models = {'greedy', 'softmax', 'rs'};
% labels = {'Greedy', 'Softmax', 'RS'};
% hexColors = {'#56B4E9', '#009E73', '#065B8D'};
% colors = cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x',[1 3])/255, hexColors, 'UniformOutput', false);
% 
% regionNames = fieldnames(brain_data);
% nRegions = numel(regionNames);
% 
% % --- Layout setup ---
% nCols = ceil(sqrt(nRegions));
% nRows = ceil(nRegions / nCols);
% 
% figure('Color','w','Position',[100 100 1600 900]);
% t = tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
% 
% for i = 1:nRegions
%     nexttile;
%     region = regionNames{i};
%     data = cell(1, numel(models));
% 
%     % collect non-NaN values per model
%     for m = 1:numel(models)
%         mdl = models{m};
%         if isfield(brain_data.(region), mdl)
%             vals = brain_data.(region).(mdl);
%             vals = vals(~isnan(vals));
%             data{m} = vals(:);
%         else
%             data{m} = [];
%         end
%     end
% 
%     % skip empty region
%     if all(cellfun(@isempty, data)), continue; end
% 
%     % --- boxplot ---
%     boxplot(cell2mat(data), repelem(1:3, cellfun(@numel, data)), ...
%         'Labels', labels, 'Colors', 'k', 'Symbol', '', 'BoxStyle', 'outline');
%     hold on;
% 
%     % --- scatter overlay ---
%     for m = 1:numel(models)
%         x = repmat(m, numel(data{m}), 1) + (rand(size(data{m})) - 0.5)*0.2;
%         scatter(x, data{m}, 20, 'filled', ...
%             'MarkerFaceColor', colors{m}, 'MarkerFaceAlpha', 0.7);
%     end
% 
%     % --- aesthetics ---
%     title(strrep(region, '_', ' '), 'Interpreter', 'none', 'FontSize', 10);
%     ylabel('R^2');
%     ylim([0 0.07]);
% 
%     ax = gca;
%     box off;
%     ax.TickDir = 'out';
%     ax.LineWidth = 1;
% end
% 
% title(t, 'R^2 per brain area', ...
%     'FontSize', 14, 'FontWeight', 'bold');
% 
% exportgraphics(gcf, fullfile(main_folder, 'R2_each_area_across_patients.pdf'), ...
%     'ContentType', 'vector', 'Resolution', 300);


%%
models = {'greedy', 'softmax', 'rs'};
labels = {'greedy', 'softmax', 'RS'};
hexColors = {'#56B4E9', '#009E73', '#065B8D'};
colors = cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x',[1 3])/255, hexColors, 'UniformOutput', false);

allR2 = struct();
allLocs = struct();

% --- NEW: collect total channels per area across patients (sig + non-sig)
allTotalLocs = {};

ptFields = fieldnames(database);
for i = 1:numel(ptFields)
    ptName = ptFields{i};
    data = database.(ptName);

    if ~isfield(data, 'results')
        warning('Skipping %s: missing results', ptName);
        continue;
    end
    results = data.results;

    % Make sure anatomical labels exist
    if ~isfield(data, 'NewAnatomicalLocs')
        warning('Skipping %s: missing NewAnatomicalLocs', ptName);
        continue;
    end
    anatomicalLocs = data.NewAnatomicalLocs;

    % --- NEW: append ALL channels' anatomical labels for totals
    if isstring(anatomicalLocs), anatomicalLocs = cellstr(anatomicalLocs); end
    allTotalLocs = [allTotalLocs; anatomicalLocs(:)];

    for m = 1:numel(models)
        mdl = models{m};
        if ~isfield(results, mdl) || ...
           ~isfield(results.(mdl), 'cvPseudoR2') || ...
           ~isfield(results.(mdl), 'pvals')
            continue;
        end

        R2 = results.(mdl).cvPseudoR2;
        p  = results.(mdl).pvals;

        sigIdx = find(p < 0.05 & ~isnan(R2));
        sigR2  = R2(sigIdx);
        sigLoc = anatomicalLocs(sigIdx);

        if isempty(sigR2), continue; end

        if ~isfield(allR2, mdl)
            allR2.(mdl) = [];
            allLocs.(mdl) = {};
        end

        % Append values and corresponding locations
        allR2.(mdl)   = [allR2.(mdl); sigR2(:)];
        allLocs.(mdl) = [allLocs.(mdl); sigLoc(:)];
    end
end


% second plot

data = {allR2.greedy, allR2.softmax, allR2.rs};
locs = {allLocs.greedy, allLocs.softmax, allLocs.rs};


hexColors = {'#56B4E9', '#009E73', '#065B8D'};
colors = cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x',[1 3])/255, ...
                 hexColors, 'UniformOutput', false);

% Remove empty entries
nonEmptyIdx = ~cellfun(@isempty, data);
data = data(nonEmptyIdx);
locs = locs(nonEmptyIdx);
labels = labels(nonEmptyIdx);
colors = colors(nonEmptyIdx);

% --- Prepare grouping for boxplot ---
groupVec = [];
dataVec  = [];
for m = 1:numel(data)
    dataVec  = [dataVec; data{m}(:)];
    groupVec = [groupVec; repmat(m, numel(data{m}), 1)];
end

% --- Figure layout ---
figure('Color','w','Position',[300 300 800 500]);
t = tiledlayout(1,2,'TileSpacing','tight','Padding','compact');

nexttile(1); hold on
boxplot(dataVec, groupVec, ...
    'Labels', labels, 'Colors', 'k', 'Symbol', '', 'BoxStyle', 'outline');

for m = 1:numel(data)
    x = repmat(m, numel(data{m}), 1) + (rand(size(data{m})) - 0.5)*0.2;
    scatter(x, data{m}, 15, ...
        'filled', ...
        'MarkerFaceColor', colors{m}, ...
        'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', 'none');
end

ylabel('R^2');
title('models R^2');
ylim([0 0.15]);
axis square;

ax = gca;
box off;
ax.TickDir = 'out';
ax.LineWidth = 1;

% --- Kruskal-Wallis ---
[pKW,~,stats] = kruskalwallis(dataVec, groupVec, 'off');
disp(['Kruskal-Wallis p = ' num2str(pKW)]);

% --- Post-hoc Wilcoxon pairwise comparisons ---
pairs = nchoosek(1:numel(data), 2);
yMax = max(dataVec) * 1.1;
yStep = range(dataVec) * 0.15;
for i = 1:size(pairs,1)
    x = data{pairs(i,1)};
    y = data{pairs(i,2)};
    p = ranksum(x, y);
    disp([labels{pairs(i,1)} ' vs ' labels{pairs(i,2)} ' p = ' num2str(p)]);
    
    a = pairs(i,1);
    b = pairs(i,2);
    yBar = yMax + (i-1)*yStep;
    lineHeight = yStep * 0.2;
    
    plot([a b], [yBar yBar], 'k', 'LineWidth', 0.5);
    plot([a a], [yBar - lineHeight, yBar], 'k', 'LineWidth', 0.5);
    plot([b b], [yBar - lineHeight, yBar], 'k', 'LineWidth', 0.5);

    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
    
    text(mean([a b]), yBar + yStep/6, stars, ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

nexttile(2); hold on

% --- Combine Left/Right labels for SIG channels ---
mergedLocs = cell(size(locs));
for m = 1:numel(locs)
    tmp = locs{m};
    tmp = regexprep(tmp, '^(Left|Right|left|right)\s+', '');
    mergedLocs{m} = tmp;
end

% --- Count occurrences of each merged anatomical label per model (SIG only) ---
allAreas = unique(cat(1, mergedLocs{:}));
areaCounts = zeros(numel(allAreas), numel(mergedLocs));
for m = 1:numel(mergedLocs)
    [u,~,idx] = unique(mergedLocs{m});
    counts = accumarray(idx, 1);
    for j = 1:numel(u)
        matchIdx = strcmp(allAreas, u{j});
        areaCounts(matchIdx,m) = counts(j);
    end
end

% keep only areas present in ALL models (as you specified)
keepIdx = all(areaCounts >= 1, 2);
allAreas = allAreas(keepIdx);
areaCounts = areaCounts(keepIdx, :);

% --- Sort by total count (descending) BEFORE normalization (unchanged) ---
[~, sortIdx] = sort(sum(areaCounts,2), 'descend');
allAreas = allAreas(sortIdx);
areaCounts = areaCounts(sortIdx,:);

% --- Build total channel counts per (merged) area across patients ---
mergedTotalLocs = regexprep(allTotalLocs, '^(Left|Right|left|right)\s+', '');
[uTot,~,idxTot] = unique(mergedTotalLocs);
countsTot = accumarray(idxTot, 1);

totalCounts = zeros(numel(allAreas),1);
for k = 1:numel(allAreas)
    hit = strcmp(uTot, allAreas{k});
    if any(hit)
        totalCounts(k) = countsTot(hit);
    else
        totalCounts(k) = 0; % safety
    end
end

% --- Normalize per area by total channels (sig + non-sig) ---
normAreaCounts = bsxfun(@rdivide, areaCounts, totalCounts);

% --- Clean sort: order areas by the max of each group's bars (desc) ---
[~, sortIdx2] = sort(max(normAreaCounts, [], 2), 'descend');
normAreaCounts = normAreaCounts(sortIdx2, :);
allAreas       = allAreas(sortIdx2);

% --- Plot grouped bars (now normalized) ---
b = bar(normAreaCounts, 'grouped', 'BarWidth', 0.9);
for m = 1:numel(mergedLocs)
    b(m).FaceColor = colors{m};
    b(m).EdgeColor = 'none';
end

% --- Axis + style ---
set(gca, 'XTick', 1:numel(allAreas), ...
    'XTickLabel', allAreas, ...
    'XTickLabelRotation', 45, ...
    'FontSize', 8, ...
    'TickDir', 'out', 'LineWidth', 1);


ylabel('proportion of channels (sig / total) (a.u.)');
title('brain areas with significant channels');

legend(labels, 'Location', 'northeastoutside');
axis square;
box off;

% --- Save figure (unchanged name/path) ---
exportgraphics(gcf, fullfile(main_folder, ...
    'R2_all_channels_sigAcrossPatients.pdf'), ...
    'ContentType', 'vector', 'Resolution', 300);

