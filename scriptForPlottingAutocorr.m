function [f,S]=scriptForPlottingAutocorr(psth,autocorr_out,usel,uses,i,lagsWindow)

% i=2;
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% usel=[freqs+0.05];
% usel=[0.05];
% usel=[0.05];
% uses=[12];

l=psth.unitLED{1};
s=psth.unitStimcond{1};
useTrials1=[];
for j=1:length(usel)
    useTrials1=[useTrials1 find(l>=usel(j)-0.001 & l<=usel(j)+0.001)];
end
useTrials1=sort(useTrials1);
useTrials2=[];
for j=1:length(uses)
    useTrials2=[useTrials2 find(s>=uses(j)-0.001 & s<=uses(j)+0.001)];
end
useTrials2=sort(useTrials2);
useTrials=useTrials1(ismember(useTrials1,useTrials2));

% figure(); 
temp=psth.psths{i};
% plot(downSampAv(psth.t,3),downSampAv(nanmean(temp(useTrials,:),1),3));
outS=nanmean(nanmean(temp(useTrials,psth.t>=0.25 & psth.t<=1.75),1),2);

% figure(); 
temp=autocorr_out.autocorrs{i};
a=nanmean(temp(useTrials,:),1);
% plot(autocorr_out.lags.*20,a,'Color','k');
% ylim([0 max(a(autocorr_out.lags.*20>0))]);

params.Fs=50;
params.tapers=[5 9];
params.trialave=1;
params.fpass=[2.5 20];
temp(temp<10^-5)=0;
% [S,f]=mtspectrumpb(nanmean(temp(useTrials,autocorr_out.lags.*20>lagsWindow(1) & autocorr_out.lags.*20<lagsWindow(2)),1),params);
% if any(temp(useTrials,autocorr_out.lags.*20>lagsWindow(1) & autocorr_out.lags.*20<lagsWindow(2))>0)
    [S,f]=mtspectrumpb(temp(useTrials,autocorr_out.lags.*20>lagsWindow(1) & autocorr_out.lags.*20<lagsWindow(2))',params);
% else
%     f=[params.fpass(1):0.5:params.fpass(2)];
%     S=zeros(size(f));
% end
% figure(); plot(f,S,'g');
% disp('hi')