function [x,psth_mu,avSpec,f,ps]=gammaInThalamus_fromPSTH(spikes,allAssigns)

bin=1;
% trialDuration=3.5;
% trialDuration=4;
trialDuration=5;
% window=[0 4]; % in s from trial osnet
% window=[0.8 3.1]; % in s from trial osnet
% window=[0.3 1.8]; % in s from trial osnet % Best for m206
% window=[0 3.5]; % in s from trial osnet
% window=[0.8 3.5]; % in s from trial osnet % m116
% window=[1.2 3.2]; % in s from trial osnet % m88A
window=[0.5 4]; % in s from trial osnet % m88A, m90b

% Make specgram from MU
l=unique(spikes.trials);
[~,~,~,x,~,psth_mu]=psth_wStd_trialByTrial(spikes,bin,0,trialDuration,length(l),l);
disp('Done with trial-by-trial PSTH');
% params.tapers=[0.9 trialDuration 0];
% params.tapers=[1 3];
params.tapers=[5 6];
% params.tapers=[2 8];
params.Fs=1/(bin/1000);
params.trialave=0;
% [S,t,f]=mtspecgrampb(psth_mu',[1 0.05],params);
[S,t,f]=mtspecgrampb(psth_mu',[0.5 0.05],params);

ps=nan(size(S,3),length(f));
for i=1:size(S,3)
    ps(i,:)=nanmean(S(t>=window(1) & t<=window(2),:,i),1);
end

concS=nan(size(S,1)*size(S,3),size(S,2));
k=1;
for i=1:size(S,3)
    concS(k:k+size(S,1)-1,:)=S(:,:,i);
    k=k+size(S,1);
end

fax=1:length(f);
subf=f(fax);
flabel=cell(1,length(fax));
for i=1:length(fax)
    flabel{i}=num2str(subf(end-(i-1)));
end

figure(); 
procS=flipud(S(:,:,1)');
imagesc(t,fax,procS);
set(gca,'Ytick',fax);
set(gca,'YtickLabel',flabel);
% ylim([fax(end-100) fax(end)]);
% ylim([fax(end-70) fax(end-50)]);
ylim([fax(1) fax(end-50)]);

long_t=linspace(0,max(t)*size(S,3),size(concS,1));
figure(); 
procS=flipud(concS');
imagesc(long_t,fax,procS);
set(gca,'Ytick',fax);
set(gca,'YtickLabel',flabel);
% ylim([fax(end-100) fax(end)]);
% ylim([fax(end-70) fax(end-50)]);
ylim([fax(1) fax(end-50)]);

figure(); 
procS=flipud(concS(:,f>=50 & f<=70)');
% procS=flipud(concS(:,f>=40 & f<=80)');
imagesc(long_t,f(f>=50 & f<=70),procS);
avSpec=nanmean(concS,1);

figure(); 
semilogy(f,avSpec);

