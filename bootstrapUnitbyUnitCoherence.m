function [C,f,phi]=bootstrapUnitbyUnitCoherence(dLGNunitpsths,V1unitpsths)

binsize=1; % in ms
params.pad=0;
% params.tapers=[0.2857 10.5 0];
% params.tapers=[0.2857 10.5 0];
params.tapers=[15 18];
params.fpass=[1 70];
params.trialave=1;
params.Fs=1000/binsize;

takeFraction=0.95;
nTrials=100;
avdLGNunitresponse=zeros(length(dLGNunitpsths),size(dLGNunitpsths{1},2));
for j=1:length(dLGNunitpsths)
    avdLGNunitresponse(j,:)=nanmean(dLGNunitpsths{j},1);
end
medLGN=median(avdLGNunitresponse,1);
avV1unitresponse=zeros(length(V1unitpsths),size(V1unitpsths{1},2));
for j=1:length(V1unitpsths)
    avV1unitresponse(j,:)=nanmean(V1unitpsths{j},1);
end
LGNforCoh=repmat(medLGN,size(avV1unitresponse,1),1);
    
%     bootdLGN=zeros(nTrials,length(nanmean(dLGN,1)));
%     bootV1=zeros(nTrials,length(nanmean(V1,1)));
%     for i=1:nTrials
%         r=randi([1 size(dLGN,1)],1,floor(size(dLGN,1)*takeFraction));
%         bootdLGN(i,:)=nanmean(dLGN(r,:),1);
%         bootV1(i,:)=nanmean(V1(r,:),1);
%     end

[C,phi,~,~,~,f]=coherencypb(LGNforCoh(:,500:9500)',avV1unitresponse(:,500:9500)',params,0);
figure();
semilogx(f,C);