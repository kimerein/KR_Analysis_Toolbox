function plot_autocorr(spikes,useAssign,ledValue)

stimWindow=[1 3];
corr_bin_size=2; % in ms;
maxlag=2; % in s;

if ~isempty(useAssign)
    spikes=filtspikes(spikes,0,'assigns',useAssign);
end
if ~isempty(ledValue)
    temp=[];
    temp1=[];
    for i=1:length(ledValue)
        spikes=makeTempField(spikes,'led',ledValue(i));
        temp(i,:)=spikes.temp;
        temp1(i,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    spikes=filtspikes(spikes,0,'temp',1);
end
s=spikes.spiketimes;

spiketimes=s(s>=stimWindow(1) & s<=stimWindow(2));
[cross,lags]=pxcorr(spiketimes,spiketimes,round(1000/corr_bin_size),maxlag);
cross(find(lags==0))=0;
ymax=max(cross)+1;
% ymin=min(cross(lags<=-10 & lags>=-20));
f1=figure();
hold on;
bb=bar(lags*1000,cross,1.0);
evoked.x=lags*1000;
evoked.y=cross;
set(bb,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
set(gca,'XLim',maxlag*1000*[-0.05 0.05]);
set(gca,'YLim',[0 ymax]);
xlabel('Time lag (ms)');
ylabel('Autocorrelation (Hz)');
title('Evoked');