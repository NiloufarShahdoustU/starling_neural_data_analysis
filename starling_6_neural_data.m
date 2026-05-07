% in this code I'm trying to read some patients data

% trialStart: 1 , by fixation point: 500ms
% cardShow: 2 , 1000ms
% instructionMessage: 3 , no delay
% flipSpace: 4 , no delay <3000ms
% choiceAndFeedback: 5 , feedback is shown for 2000ms
% totalReward: 6 , shown for 1000ms

clc;
clear;
close all;
% add new patients here
differentPatients = {'202514', '202518', '202521', '202522' , '202601'};

%% all these times are in ms
cueStart = 100;
cueEnd = 900;

flipStart = 500;
flipEnd = 500;

choiceStart = 800;
choiceEnd = 200;

totalRewardStart = 200;
totalRewardEnd = 800;

%% Main folder
input_folder = fullfile('\\155.100.91.44\d\Data\Nill\starling\raw');

d = dir(input_folder);

isub = [d(:).isdir]; 
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = [];
ptIDs = string(subFolders);

%%

% ATTENTION: when reading ns2 for the new system (from patient 202514 on)
% when we open ns2, and Data, there are 4 cells. we need to concatenate
% those data from those cells to have a proper complete data. 

% for p = 1:numel(ptIDs)

for p = 9:9

    ptID = ptIDs{p};
    fprintf('\n--- Processing ptID: %s ---\n', ptID);
    
    input_folder_pt = fullfile(input_folder, ptID); 

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
    selectedChans = find(~isECoG);
    selectedChans = selectedChans(1:end-1); 
    SelectedAnatomicalLoc = anatomicalLocs(selectedChans);

    nChans = length(selectedChans);
    LFPData = [];
    % notch filter on 60 Hz
    [b1,a1] = iirnotch(60/(original_freq/2),(60/(original_freq/2))/25);
    
    for ch = 1:nChans
        % weird 202514 and 202521 case:
        if ismember(ptID, {'202514','202521'})
            % concatenate all cells in NS2.Data along time dimension
            tmpData = [];
            for c = 1:numel(NS2.Data)
                tmpData = [tmpData, double(NS2.Data{c}(selectedChans(ch),:))];
            end
            tmp = filtfilt(b1, a1, tmpData);
            LFPData(ch,:) = tmp;
            clear tmpData
        % other cases:
        else
            tmp = filtfilt(b1, a1, double(NS2.Data(selectedChans(ch),:)));      
            LFPData(ch,:) = tmp;
        end  
        clear tmp
    end
   

   % reading task data:
    bhvFiles = dir(fullfile(input_folder_pt, 'task_data*.csv'));
    bhvFile = fullfile(bhvFiles(1).folder, bhvFiles(1).name);
    bhvData = readtable(bhvFile);

    % reading eventTimes:
    eventTimesFile = fullfile(input_folder_pt, 'eventTimes.mat');
    eventTimes = load(eventTimesFile);


end




%% DEBUG
 


