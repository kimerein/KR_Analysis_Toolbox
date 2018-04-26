function [C,f,phi]=bootstrapCoherence(dLGN,V1)

binsize=1; % in ms
differentTrialSets=0;
params.pad=0;
% params.tapers=[0.2857 10.5 0];
% params.tapers=[0.2857 10.5 0];
params.tapers=[15 18];
params.fpass=[1 70];
params.trialave=1;
params.Fs=1000/binsize;
takeWindow=[0 10.5];

takeFraction=0.95;
nTrials=100;
bootdLGN=zeros(nTrials,length(nanmean(dLGN,1)));
bootV1=zeros(nTrials,length(nanmean(V1,1)));
for i=1:nTrials
    disp(i);
    if differentTrialSets==0
        r=randi([1 size(dLGN,1)],1,floor(size(dLGN,1)*takeFraction));
        bootdLGN(i,:)=nanmean(dLGN(r,:),1);
        bootV1(i,:)=nanmean(V1(r,:),1);
    else
        r1=randi([1 size(dLGN,1)],1,floor(size(dLGN,1)*takeFraction));
        r2=randi([1 size(V1,1)],1,floor(size(V1,1)*takeFraction));
        bootdLGN(i,:)=nanmean(dLGN(r1,:),1);
        bootV1(i,:)=nanmean(V1(r2,:),1);
    end
end

[C,phi,~,~,~,f]=coherencypb(bootdLGN(:,takeWindow(1)*1000:takeWindow(2)*1000)',bootV1(:,takeWindow(1)*1000:takeWindow(2)*1000)',params,0);
figure();
semilogx(f,C);