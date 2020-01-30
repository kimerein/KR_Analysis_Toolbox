function compareVisStimWithResponseCoherence(noTheta_psth,dLGNpsth,noThetaTrialVal,noThetaTrials,stimFreq,Fs)

ds=1;

t=noTheta_psth.t;

% set up vis stim
y=sin(2*pi*stimFreq.*t(t>0 & t<=2));

trialsTogether=[];
temp=noTheta_psth.psths{1};
l=noTheta_psth.unitLED{1};
trialSum=zeros(size(downSampMatrix(temp(ismember(l,[stimFreq]),:),ds)));
for i=1:length(noTheta_psth.psths)
    temp=noTheta_psth.psths{i};
    l=noTheta_psth.unitLED{i};
    outs(i,:)=nanmean(temp(ismember(l,[stimFreq]),:),1);
    trialSum=trialSum+downSampMatrix(temp(ismember(l,[stimFreq]),:),ds);
    trialsTogether=[trialsTogether; temp(ismember(l,[stimFreq]),:)];
end
% figure();
% for i=1:size(trialSum,1)
%     plot(downSampAv(noTheta_psth.t,ds),(i-1)*100+trialSum(i,:));
%     hold all;
% end

% return

takeunitshere=[5 11];
temp=dLGNpsth.psths{1};
l=dLGNpsth.unitLED{1};
trialSum2=zeros(size(downSampMatrix(temp(ismember(l,[stimFreq]) & noThetaTrials'==noThetaTrialVal,:),ds)));
for i=1:length(takeunitshere)
    temp=dLGNpsth.psths{i};
    l=dLGNpsth.unitLED{i};
    out2(i,:)=nanmean(temp(ismember(l,[stimFreq]) & noThetaTrials'==noThetaTrialVal,:),1);
    trialSum2=trialSum2+downSampMatrix(temp(ismember(l,[stimFreq]) & noThetaTrials'==noThetaTrialVal,:),ds);
    trialsTogether=[trialsTogether; temp(ismember(l,[stimFreq]) & noThetaTrials'==noThetaTrialVal,:)];
end

tog=[outs; out2];

s=15;
for i=1:size(tog,1)
    tog(i,:)=smooth(tog(i,:),s);
end

figure(); plot(downSampAv(dLGNpsth.t,ds),smooth(nanmean(downSampMatrix(tog,ds),1),1)); title('Average PSTH');

% figure(); plot(downSampAv(dLGNpsth.t,ds),trialSum'); title('Single trials summed PSTH');
% 
% figure(); plot(downSampAv(dLGNpsth.t,ds),trialSum2'); title('Single trials summed PSTH');

return

% params.tapers=[3 5];
params.tapers=[1 3];
% params.tapers=[18 21];
params.Fs=Fs;
params.fpass=[1 50];
params.pad=0;
params.trialave=1;

coh=[];
ph=[];
tdelay=[];
% [C,phi,S12,S1,S2,f]=coherencypb(nanmean(repmat(y,size(trialsTogether,1),1),1)',nanmean(trialsTogether(:,t>1 & t<=3),1)',params);
[C,phi,S12,S1,S2,f]=coherencypb(repmat(y,size(trialsTogether,1),1)',trialsTogether(:,t>1 & t<=3)',params);

figure();
plot(f,nanmean(C,2));
title('Coherence');

figure();
plot(f,nanmean(phi,2));
% hold on;
% plot(f,nanmean(phi,2)+nanstd(phi,[],2));
title('Phase');