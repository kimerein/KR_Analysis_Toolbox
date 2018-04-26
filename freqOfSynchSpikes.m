function [take,allTrials_difference,takeSpecgram]=freqOfSynchSpikes(synchPSTH,timewindow,nonprefSynchPSTH)

if iscell(synchPSTH)
    synchPSTH=synchPSTH{1};
    nonprefSynchPSTH=nonprefSynchPSTH{1};
end

trialDuration=14.5;
stimWindow=[4 6.5];
ds=1;
times=linspace(0,trialDuration,size(synchPSTH,2));

params.tapers=[3 10];
params.Fs=1./(trialDuration/size(synchPSTH,2));
params.fpass=[1 20];
params.trialave=0;
movingwin=[1 0.05];
% movingwin=[1 0.05];

[S,f]=mtspectrumpb(synchPSTH(:,times>=timewindow(1) & times<=timewindow(2))',params);
[Sev,fev]=mtspectrumpb(synchPSTH(:,times>=stimWindow(1) & times<=stimWindow(2))',params);
[Snonpref,fnonpref]=mtspectrumpb(nonprefSynchPSTH(:,times>=timewindow(1) & times<=timewindow(2))',params);

params.tapers=[3 10];
params.Fs=1./(trialDuration/size(synchPSTH,2));
params.fpass=[1 20];
params.trialave=1;
[S_pref,t_pref,f_pref]=mtspecgrampb(synchPSTH',movingwin,params);
[S_nonpref,t_nonpref,f_nonpref]=mtspecgrampb(nonprefSynchPSTH',movingwin,params);

figure(); 
imagesc(t_pref,f_pref,S_pref');
figure(); 
imagesc(t_nonpref,f_nonpref,S_nonpref');
figure(); 
imagesc(t_nonpref,f_nonpref,S_pref'-S_nonpref');
allTrials_difference=S_pref'-S_nonpref';

% figure(); 
% scatter(nanmean(S(f>=4 & f<=8,:),1),nanmean(S(f>8 & f<=16,:),1));

LFalpha=nanmean(S(f>=4 & f<=6,:),1);
HFalpha=nanmean(S(f>=11.5 & f<=12.5,:),1);
% LFalpha=nanmean(Snonpref(fnonpref>=4 & fnonpref<=6,:),1);
% HFalpha=nanmean(Snonpref(fnonpref>=11.5 & fnonpref<=12.5,:),1);
ratio=LFalpha./HFalpha;

take=ratio<1;
% take=HFalpha<mean(HFalpha);
% take=nanmean(S(f>=2.5 & f<=3.5,:),1)>=median(nanmean(S(f>=2.5 & f<=3.5,:),1));
% take=LFalpha>mean(LFalpha);
% take=nanmean(Sev(fev>=2.5 & fev<=3.5,:),1)>=median(nanmean(Sev(fev>=2.5 & fev<=3.5,:),1));



params.trialave=1;
[S_pref,t_pref,f_pref]=mtspecgrampb(synchPSTH(take,:)',movingwin,params);
[S_nonpref,t_nonpref,f_nonpref]=mtspecgrampb(nonprefSynchPSTH(take,:)',movingwin,params);

figure(); 
imagesc(t_pref,f_pref,S_pref');
figure(); 
imagesc(t_nonpref,f_nonpref,S_nonpref');
figure(); 
imagesc(t_nonpref,f_nonpref,S_pref'-S_nonpref');
takeSpecgram=S_pref'-S_nonpref';


figure(); 
plot(downSampAv(times,ds),downSampAv(nanmean(synchPSTH(take,:),1),ds),'Color','r');

figure(); 
plot(downSampAv(times,ds),downSampAv(nanmean(nonprefSynchPSTH(take,:),1),ds),'Color','b');

figure(); 
plot(downSampAv(times,ds),downSampAv(nanmean(synchPSTH(take,:),1),ds),'Color','r');
hold on;
plot(downSampAv(times,ds),downSampAv(nanmean(nonprefSynchPSTH(take,:),1),ds),'Color','b');
% figure();
% scatter(nanmean(S(f>=4 & f<=6,:),1),nanmean(Sev(fev>=2.5 & fev<=3.5,:),1));

figure();
scatter(ratio,nanmean(Sev(fev>=2.5 & fev<=3.5,:),1));
disp('done');

figure(); 
b=bandPassLFP(synchPSTH(take,:),1/(times(2)-times(1)),11.5,12.5,0);
runningy=0;
for i=1:size(b,1)
    plot(times,b(i,:)+runningy);
    hold all;
    runningy=runningy+max(b(i,:));
end