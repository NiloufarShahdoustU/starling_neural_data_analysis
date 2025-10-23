clc;
clear;
close all;


nevList = dir('\\155.100.91.44\d\Data\Nill\starling\raw\202514\*.nev');



if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end

%% load neural data from ns5 to get photodiode
[nevPath,nevName,nevExt] = fileparts(nevFile);
NS5 = openNSx(fullfile(nevPath,[nevName '.ns5']));
%% photodiode
close all;
plot(NS5.Data(2, :))
