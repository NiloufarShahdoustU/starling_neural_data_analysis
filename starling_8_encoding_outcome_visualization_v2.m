clc;
clear;
close all

main_folder = '\\155.100.91.44\d\Code\Nill\Starling_neural_data\starling_8_encoding_outcome_v2';

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
models = {'greedy', 'softmax', 'rs'};
labels = {'Greedy', 'Softmax', 'RS'};
hexColors = {'#56B4E9', '#009E73', '#065B8D'};
colors = cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x',[1 3])/255, hexColors, 'UniformOutput', false);

regionNames = fieldnames(brain_data);
nRegions = numel(regionNames);

% --- Layout setup ---
nCols = ceil(sqrt(nRegions));
nRows = ceil(nRegions / nCols);

figure('Color','w','Position',[100 100 1600 900]);
t = tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');

for i = 1:nRegions
    nexttile;
    region = regionNames{i};
    data = cell(1, numel(models));
    
    % collect non-NaN values per model
    for m = 1:numel(models)
        mdl = models{m};
        if isfield(brain_data.(region), mdl)
            vals = brain_data.(region).(mdl);
            vals = vals(~isnan(vals));
            data{m} = vals(:);
        else
            data{m} = [];
        end
    end
    
    % skip empty region
    if all(cellfun(@isempty, data)), continue; end
    
    % --- boxplot ---
    boxplot(cell2mat(data), repelem(1:3, cellfun(@numel, data)), ...
        'Labels', labels, 'Colors', 'k', 'Symbol', '', 'BoxStyle', 'outline');
    hold on;
    
    % --- scatter overlay ---
    for m = 1:numel(models)
        x = repmat(m, numel(data{m}), 1) + (rand(size(data{m})) - 0.5)*0.2;
        scatter(x, data{m}, 20, 'filled', ...
            'MarkerFaceColor', colors{m}, 'MarkerFaceAlpha', 0.7);
    end
    
    % --- aesthetics ---
    title(strrep(region, '_', ' '), 'Interpreter', 'none', 'FontSize', 10);
    ylabel('R^2');
    ylim([0 0.07]);
    
    ax = gca;
    box off;
    ax.TickDir = 'out';
    ax.LineWidth = 1;
end

title(t, 'R^2 per brain area', ...
    'FontSize', 14, 'FontWeight', 'bold');

exportgraphics(gcf, 'R2_each_area_across_patients.pdf', 'ContentType', 'vector', 'Resolution', 300);

