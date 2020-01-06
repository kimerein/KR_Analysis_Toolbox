function doAnalysisFromSpikesToF1(expt,spikes,noThetaTrials,take,outputDir,tookFileInds)

% tookFileInds=[1:155];
trialDuration=14.5;
% outputDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Ntsr1 opto stim\Mawake397\dLGN T1 LED sines';
ledOffVal=[0 0.05];
% take=[];
putStimAt=4; % in seconds

detectfilenames=expt.files.names(tookFileInds);
spikes=fixSpikesAfterO2Sort(spikes,expt,detectfilenames);
save([outputDir '\' spikes.info.spikesfile '.mat'],'spikes');
temp=spikes.labels(:,1);
useAssigns=temp(spikes.labels(:,2)==2); % 2 means good unit
disp(useAssigns);
psth=measureAllUnitsResponseProperties(spikes,useAssigns,[0 trialDuration]);
save([outputDir '\psth_notAlignedToStim.mat'],'psth');
psth=alignToStimFromExpt(psth,expt,tookFileInds,trialDuration,putStimAt,isnan(spikes.sweeps.led));
save([outputDir '\psth_alignedToStim.mat'],'psth');
s=unique(spikes.stimcond(~isnan(spikes.stimcond)));
disp(s);
clear uses; clear uses_tri;
if isempty(take)
    a=unique(spikes.trials(~isnan(spikes.stimcond))); % use all trials
else
    a=take; % use only these trials
end
for i=1:length(s)
    uses{i}={[s(i)]};
    uses_tri{i}={a};
end
led=unique(psth.unitLED{1});
ledOnVal=led(led>=0.01);
if length(noThetaTrials)~=size(psth.psths{1},1)
    noThetaTrials=noThetaTrials(~isnan(spikes.sweeps.led));
end
doF1analysis([],[],outputDir,[],psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials);
