function params=getLFPAnalysisParams(expt)
% Extract experiment-specific parameters from experiment struct
% I believe this function has been replaced by getAnalysisParams.m in
% automatic experiment analysis package

nAnalysisBlocks=1;
filesInBlocks=[1 1];
currStimParams=expt.stimulus(1).params;
for i=2:length(expt.stimulus)
    if isequal(expt.stimulus(i).params,currStimParams)
        filesInBlocks(nAnalysisBlocks,2)=i;
    else
        nAnalysisBlocks=nAnalysisBlocks+1;
        filesInBlocks=[filesInBlocks; i i];
        currStimParams=expt.stimulus(i).params;
    end
end
stimBlocks=cell(nAnalysisBlocks,1);
for i=1:size(filesInBlocks,1)
    stimBlocks{i}=expt.stimulus(filesInBlocks(i,1)).params.StimulusName;
end
disp('These are the different stimulus parameter blocks for this experiment.');
disp(stimBlocks);
disp('Stimulus parameter blocks correspond to these DAQ files. (Col. 1 through Col. 2 for each row)');
disp(filesInBlocks);
c=input('Continue?');
params=[];
if isempty(c)
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
        s.Var2_name=expt.stimulus(filesInBlocks(i,1)).varparam(2).Name;
        s.Var2_values=expt.stimulus(filesInBlocks(i,1)).varparam(2).Values;
        params=[params; s];
    end 
else
    params=[];
end
    