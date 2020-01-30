function doAnalysisFromSpikesToF1(expt,spikes,noThetaTrials,take,outputDir,tookFileInds,onlyUseAssigns)

% tookFileInds=[1:155];
% trialDuration=14.5;
trialDuration=4;
% outputDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Ntsr1 opto stim\Mawake397\dLGN T1 LED sines';
% ledOffVal=[0 0.05];
ledOffVal=[0];
% take=[];
% putStimAt=4; % in seconds
putStimAt=1; % in seconds

detectfilenames=expt.files.names(tookFileInds);
if ~isfield(spikes,'stimcond')
    spikes=fixSpikesAfterO2Sort(spikes,expt,detectfilenames);
    save([outputDir '\' spikes.info.spikesfile '.mat'],'spikes');
end
if ~isempty(onlyUseAssigns)
    spikes=filtspikes(spikes,0,'assigns',onlyUseAssigns);
    useAssigns=unique(spikes.assigns);
else
    temp=spikes.labels(:,1);
    useAssigns=temp(spikes.labels(:,2)==2); % 2 means good unit
end
disp(useAssigns);
psth=measureAllUnitsResponseProperties(spikes,useAssigns,[0 trialDuration]);
save([outputDir '\psth_notAlignedToStim.mat'],'psth');
% psth=alignToStimFromExpt(psth,expt,tookFileInds,trialDuration,putStimAt,isnan(spikes.sweeps.led));
% save([outputDir '\psth_alignedToStim.mat'],'psth');
s=unique(spikes.stimcond(~isnan(spikes.stimcond)));
disp(s);
clear uses; clear uses_tri;
if isempty(take)
    a=unique(spikes.trials(~isnan(spikes.stimcond))); % use all trials
else
    a=take; % use only these trials
end
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
for i=1:length(freqs)
    uses{i}={[freqs(i) freqs(i)+0.05]};
    uses_tri{i}={a};
end
% for i=1:length(s)
%     uses{i}={[s(i)]};
%     uses_tri{i}={a};
% end
% led=unique(psth.unitLED{1});
% ledOnVal=led(led>=0.01);
if length(noThetaTrials)~=size(psth.psths{1},1)
    noThetaTrials=noThetaTrials(~isnan(spikes.sweeps.led));
end
save([outputDir '\noThetaTrials.mat'],'noThetaTrials');
% doF1analysis([],[],outputDir,[],psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials);
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
doF1analysis_freqs([],[],outputDir,[],psth,single(freqs),single(freqs+0.05),[],uses,uses_tri,noThetaTrials);
