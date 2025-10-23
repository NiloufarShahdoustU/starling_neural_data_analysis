function [labels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesSTARLING(patientID)
% PTTRODESBART gets electrode information from a .csv file
%
%   [labels,isECoG,isEEG,isECG,adjacentChanMat] = ptTrodesBART(patientID)
%   gets electrode information for the patient number listed in
%   [patientID].
%
%   labels contains all electrode labels.
%
%   isECoG indexes the intracranial channels.
%
%   isEEG indexes the extracranial channels.
%
%   isECG indexes the electrocardiogram channels.
%
%   anatomicalLocs is a cell containing the NMM atlas projection for each
%   ECoG electrode. If this data is not found in the directory specified,
%   an empty cell array with length = sum(isECoG) will be output in order
%   to preserve the functionality of dependents.
%       - As of 20230228 we're replacing 'Cerebral White Matter' electrode
%       labels with the second highest label probability.
%
%   adjacentChansMat produces a connectivity matrix for intracranial sEEG
%   channels.
%       IMPORTANT NOTE:: This matrix connects adjacent channels of 10
%       channel sEEG leads, skipping every tenth electrode. i.e. We assume
%       adjacent sEEG channels are connected, but not across sEEG
%       electrodes. This matrix can be used with clust_perm2.m.
%

% author: [EHS20190125]
% updated to remove Cerebral White Matter label 20230228


if ~isstr(patientID)
    patientID = num2str(patientID);
end


% [20190315] getting anatomical labels, if they exist.
trodeLocs = fullfile('\\155.100.91.44\d','Data',['UIC' patientID],'Imaging','Registered','ChannelMap.mat');
if exist(trodeLocs,'file')
    load(trodeLocs)

    % adjusting variable names for the newer GUI
    if exist('ChanMap','var')
        fprintf('no probabilistic projections for patient %s',patientID)
        anatomicalLocs = ChanMap.ElecNMMProj;
    elseif exist('ElecNMMProj','var')
        fprintf('no probabilistic projections for patient %s',patientID)
        anatomicalLocs = ElecNMMProj;
    else
        tmpLocs = ElecAtlasProbProj(:,1);
    end

else
    fprintf('\nNo anatomical information available for patient %s',patientID)
    try
        load(fullfile('\\155.100.91.44\d','Data',['CS' patientID],'Imaging','Registered','Electrodes.mat'));
    catch
        load(fullfile('\\155.100.91.44\d','Data',['UIC' patientID],'Imaging','Registered','Electrodes.mat'));
    end

    % pre-20230228 code::
    %     anatomicalLocs = cell(size(ElecNMMProjRaw,1));
    %     anatomicalLocs(:) = {strings};

    tmpLocs = ElecAtlasProbProj(:,1);

end

%
% whichChannelMap = 2;
% switch whichChannelMap
%     case 2


%% [20230309] for looping over channels, is max the right operation, or is length? 

try
    whichChanMapStr = 'ChanMap2';
    % converting channel mappinns to a list of labels.
    if exist('ChanMap','var')
        chans = ChanMap.ChannelMap2(~isnan(ChanMap.ChannelMap2));
        NaC='NaC';
        for c = max(chans):-1:1
            % neurologist labels
            if sum(ismember(chans,c)) == 1
                labels(c,1) = ChanMap.LabelMap(ChanMap.ChannelMap2==c);
            else
                labels(c,1) = [NaC num2str (c)];
            end

            % anatomical labels
            if ~isempty(tmpLocs{c,1})
                anatomicalLocs{c} = tmpLocs{c,1}{1,1};
                anatomicalLocProbs(c) = tmpLocs{c,1}{1,2};
                if contains(anatomicalLocs(c),'Cerebral White Matter')
                    anatomicalLocs{c} = tmpLocs{c,1}{2,1};
                    anatomicalLocProbs(c) = tmpLocs{c,1}{1,2};
                end
            else
                anatomicalLocs{c} = {[NaC num2str(c)]};
                anatomicalLocProbs(c) = 0;
            end
        end
    else
        chans = ChannelMap2(~isnan(ChannelMap2));
        NaC= 'NaC';
        for c = max(chans):-1:1
            % neurologist labels
            if sum(ismember(chans,c)) == 1
                labels(c,1) = LabelMap(ChannelMap2==c);
            else
                labels(c,1) = [NaC num2str(c)];
            end

            % anatomical labels
            if ~isempty(tmpLocs{c,1})
                anatomicalLocs{c} = tmpLocs{c,1}{1,1};
                anatomicalLocProbs(c) = tmpLocs{c,1}{1,2};
                if contains(anatomicalLocs(c),'Cerebral White Matter')
                    anatomicalLocs{c} = tmpLocs{c,1}{2,1};
                    anatomicalLocProbs(c) = tmpLocs{c,1}{1,2};
                end
            else
                anatomicalLocs{c} = {[NaC num2str(c)]};
                anatomicalLocProbs(c) = 0;
            end
        end
    end

    % label labels
    isECoG = false(length(labels),1);
    isEEG = false(length(labels),1);
    isECG = ~cellfun(@isempty,strfind(labels,'ECG'));

    % labeling labels.
    for ll = 1:length(labels)
        if (labels{ll}(1)=='R' || labels{ll}(1)=='L' || strcmp(labels{ll}(1:2),'AN') || strcmp(labels{ll}(1:2),'PO') || strcmp(labels{ll}(1:2),'GR') || strcmp(labels{ll}(1),'G') || strcmp(labels{ll}(1:2),'SU') || strcmp(labels{ll}(1:2),'OF')) && ~isECG(ll)
            isECoG(ll) = true;
        elseif (~isECoG(ll) && ~isECG(ll) && ~strcmp(labels{ll},'NaC') && ~strcmp(labels{ll},'NaN'))
            isEEG(ll) = true;
        end
        % accounting for micro labels. 
        if labels{ll}(1)=='m'
            isECoG(ll) = false;
        end
    end

    % setting up adjacent channel matrix for cluster correction.
    matSize = sum(isECoG);
    adjacentChanMat = false(matSize);
    for ch = 1:matSize
        if ~isequal(mod(ch,10),0)
            adjacentChanMat(ch,ch+1) = true;
        end
    end
    % case 1
catch
    whichChanMapStr = 'ChanMap1';
    % converting channel mappinns to a list of labels.
    if exist('ChanMap','var')
        chans = ChanMap.ChannelMap1(~isnan(ChanMap.ChannelMap1));
        NaC='NaC';
        for c = max(chans):-1:1
            % neurologist labels
            if sum(ismember(chans,c)) == 1
                labels(c,1) = ChanMap.LabelMap(ChanMap.ChannelMap1==c);
            else
                labels(c,1) = {[NaC num2str(c)]};
            end

            % anatomical labels
            if ~isempty(tmpLocs{c,1})
                anatomicalLocs{c} = tmpLocs{c,1}{1,1};
                anatomicalLocProbs(c) = tmpLocs{c,1}{1,2};
                if contains(anatomicalLocs(c),'Cerebral White Matter')
                    anatomicalLocs{c} = tmpLocs{c,1}{2,1};
                    anatomicalLocProbs(c) = tmpLocs{c,1}{1,2};
                end
            else
                anatomicalLocs{c} = {[NaC num2str(c)]};
                anatomicalLocProbs(c) = 0;
            end
        end
    else
        chans = ChannelMap1(~isnan(ChannelMap1));
        NaC='NaC';
        for c = max(chans):-1:1
            % neurologist labels
            if sum(ismember(chans,c)) == 1
                labels(c,1) = LabelMap(ChannelMap1==c);
            else
                labels(c,1) = {[NaC num2str(c)]};
            end

            % anatomical labels
            if ~isempty(tmpLocs{c,1})
                anatomicalLocs{c} = tmpLocs{c,1}{1,1};
                anatomicalLocProbs(c) = tmpLocs{c,1}{1,2};
                if contains(anatomicalLocs(c),'Cerebral White Matter')
                    anatomicalLocs{c} = tmpLocs{c,1}{2,1};
                    anatomicalLocProbs(c) = tmpLocs{c,1}{1,2};
                end
            else
                anatomicalLocs{c} = {[NaC num2str(c)]};
                anatomicalLocProbs(c) = 0;
            end
        end
    end

    % label labels
    labels(cellfun(@isempty,labels)) = {['NaC' num2str(rand(1))]};
    isECoG = false(length(labels),1);
    isEEG = false(length(labels),1);
    isECG = ~cellfun(@isempty,strfind(labels,'ECG'));

    % labeling labels.
    for ll = 1:length(labels)
        if (labels{ll}(1)=='R' || labels{ll}(1)=='L' || strcmp(labels{ll}(1:2),'AN') || strcmp(labels{ll}(1:2),'PO') || strcmp(labels{ll}(1:2),'GR') || strcmp(labels{ll}(1),'G') || strcmp(labels{ll}(1:2),'SU') || strcmp(labels{ll}(1:2),'OF')) && ~isECG(ll)
            isECoG(ll) = true;
        elseif (~isECoG(ll) && ~isECG(ll) && ~strcmp(labels{ll},'NaC') && ~strcmp(labels{ll},'NaN'))
            isEEG(ll) = true;
        end
        % accounting for micro labels. 
        if labels{ll}(1)=='m'
            isECoG(ll) = false;
        end
    end

    % setting up adjacent channel matrix for cluster correction.
    matSize = sum(isECoG);
    adjacentChanMat = false(matSize);
    for ch = 1:matSize
        if ~isequal(mod(ch,10),0)
            adjacentChanMat(ch,ch+1) = true;
        end
    end

end

