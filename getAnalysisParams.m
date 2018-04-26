function params=getAnalysisParams(expt,userStimBlocks)
% Extract experiment-specific parameters from experiment struct
% and get blocks of trials with different stimulus configurations
%
% RETURNS:
% params: a structure array with an element for each stim. block
% also, fields in each element for sampling rate, stimulus type, associated
% daq file names, stimulus duration in seconds, stimulus onset in seconds
% relative to trial onset, total trial length, stimulus paramater 1 name
% and values, stimulus parameter 2 name and values

nAnalysisBlocks=1;
filesInBlocks=[1 1];
%currStimParams=expt.stimulus(1).params;
currStimParamFile=expt.stimulus(1).paramsfile;
for i=2:length(expt.stimulus)
    %if isequal(expt.stimulus(i).params,currStimParams)
    if strcmp(expt.stimulus(i).paramsfile,currStimParamFile)
        filesInBlocks(nAnalysisBlocks,2)=i;
    else
        nAnalysisBlocks=nAnalysisBlocks+1;
        filesInBlocks=[filesInBlocks; i i];
        %currStimParams=expt.stimulus(i).params;
        currStimParamFile=expt.stimulus(i).paramsfile;
    end
end
stimBlocks=cell(nAnalysisBlocks,1);
for i=1:size(filesInBlocks,1)
    stimBlocks{i}=expt.stimulus(filesInBlocks(i,1)).params.StimulusName;
end
disp('These are the different stimulus parameter blocks for this experiment.');
disp(stimBlocks);
disp('Stimulus parameter blocks correspond to these DAQ files. (Col. 1 through Col. 2 for each row)');
if ~isempty(userStimBlocks)
    disp(userStimBlocks);
else
    disp(filesInBlocks);
end
c=input('Continue? (Enter) or Specify stim. blocks manually (Type Manual)','s');
params=[];
if isempty(c)
    if ~isempty(userStimBlocks)
        filesInBlocks=userStimBlocks;
    end
else
    filesInBlocks=[];
    c=input('Stim. block. Enter daq files as, for instance, [1 2 3].','s');
    while ~isempty(c)
        filesInBlocks=[filesInBlocks; str2num(c)];
        c=input('Stim. block. Enter daq files as, for instance, [1 2 3].','s');
    end
end
for i=1:size(filesInBlocks,1)
    s=struct();
    s.samplingRate=expt.files.Fs(filesInBlocks(i,1));
    s.stimType=expt.stimulus(filesInBlocks(i,1)).params.StimulusName;
    s.daqFiles=filesInBlocks(i,1):filesInBlocks(i,2);
    s.ONlength=expt.stimulus(filesInBlocks(i,1)).params.StimDuration;
    s.ONstart=[]; % Get from photodiode data
    s.totalTrialLength=expt.files.duration(filesInBlocks(i,1));
    s.Var1_name=expt.stimulus(filesInBlocks(i,1)).varparam(1).Name;
    s.Var1_values=expt.stimulus(filesInBlocks(i,1)).varparam(1).Values;
    if isempty(s.Var1_values)
        s.Var1_values=NaN;
    end
    if length(expt.stimulus(filesInBlocks(i,1)).varparam)>1
        s.Var2_name=expt.stimulus(filesInBlocks(i,1)).varparam(2).Name;
        s.Var2_values=expt.stimulus(filesInBlocks(i,1)).varparam(2).Values;
        if isempty(s.Var2_values)
            s.Var2_values=NaN;
        end
    else
        s.Var2_name='None';
        %s.Var2_values=[];
        s.Var2_values=NaN;
    end
    params=[params; s];
end
    