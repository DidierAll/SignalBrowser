
% This script extracts time, ecg, respiration and skin conductance (EDA)
% from the Popane dataset csv files and saves them under a Reformatted
% subfodler, so that they can be read by SignalBrowser
% Original File column headers: timestamp, affect, ECG, EDA, temp, respiration, SBP, DBP, marker
% New reformatted file column headers: timestamp, ECG, respiration, EDA
% The marker column provides information about the experimental condition/stimuli which can be used to define % some Analysis Sample window in SignalBrowser 


[file_list,pathname]=uigetfile({'*.*','All Files (*.*)'},'MultiSelect', 'on');


% load metadata to find stimuli
Tstim = readtable(fullfile(pathname,'metadata.xlsx'),'Sheet','list of stimuli','Range','A2:H33');
% Add a stimuli "NaN" for "NaN"
Tstim.StimuliID(end+1) = -99;
Tstim.StimuliShortNameFile{end} = 'No Stimuli';
Tstim.StimuliValence{end} = 'NA';

pathname_new = fullfile(pathname,'Reformatted');
mkdir(pathname,'Reformatted')
if ~iscell(file_list) % in case only one file is selected
    file_list = {file_list};
end

for i = 1:length(file_list)
    %% create and save a new tabel with selected columns
    T = readtable(fullfile(pathname,file_list{i}));
    Tnew = T(:,[1 3 6 4]); % select only time, ecg, eda and respiration
    writetable(Tnew,fullfile(pathname_new,file_list{i}))

    %% create and save a table with stimilus info 
    stim = T.marker;
    stim(isnan(stim))=-99; %change NaN to -99
    ni = find(diff(stim)~=0); % find indices of change in stimuli
    ni = [ni ; length(T.timestamp)]; % add last indices of last stimuli
    i_sel = find(stim(ni)~=-99); % select those ~= -99
    if i_sel(1)>1 
        Samplet1 = T.timestamp(ni(i_sel-1)+1); % start of each stimuli 
    else % handle situation where first stimuli is not -99
        Samplet1 = T.timestamp([1;ni(i_sel(2:end)-1)+1]);
    end
    Samplet2 = T.timestamp(ni(i_sel)); % end of each stimuli
    Ns = length(i_sel);
    SampleLabel = cell(Ns,1);
    for k = 1:Ns
        SampleLabel(k) = Tstim.StimuliShortNameFile(find(Tstim.StimuliID==stim(ni(i_sel(k)))));
    end
    SampleNber = (1:Ns)';
    Nber_of_Samples = Ns;
    Tnew_Stim = table(SampleNber,Samplet1,Samplet2,SampleLabel); 
    % save stimuli info and timing in a csv file
    %file_stim = [file_list{i}(1:strfind(file_list{i},'_')) 'Sample_Parameter'];
    file_stim = [file_list{i}(1:strfind(file_list{i},'.')-1)];
    writetable(Tnew_Stim,fullfile(pathname_new,[file_stim '_StimuliInfo.csv']));
    % save stimuli info and timing as parameter mat file to be read by SignalBrowser 
    ud.Samplet1 = Samplet1; ud.Samplet2 = Samplet2; 
    ud.Nber_of_Samples = Nber_of_Samples; ud.SampleLabel = SampleLabel;
    ud.Sample_Nber = 1; 
    save(fullfile(pathname_new,[file_stim '_Stimuli_Parameter.mat']),"ud") 
end


