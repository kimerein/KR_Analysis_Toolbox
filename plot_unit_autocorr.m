function [evoked,spont,useLEDcond,maxlag,fre,avSpec,fre2,avSpec2,freSpont,avSpecSpont,freSpont2,avSpecSpont2,forfigs]=plot_unit_autocorr(spikes,useAssign,useLEDcond)

showFigs=0;
doCrossCorr=0;
showCrossCorr=0;
doFFT=1;
showPSTH=0;
corr_bin_size=3; % in ms
% maxlag=2; % in s
% maxlag=3; % in s
maxlag=3; % in s
% useLEDcond=60;
usestimcond=1:10000;
% stimWindow=[1 3];
% stimWindow=[1.2 3.2];
stimWindow=[0.5 4];
% stimWindow=[1.2 4];
% stimWindow=[0 3.5];
% stimWindow=[0 2];
% stimWindow=[0.8 3.5];
baseWindow=[0 1];
evoked=[];
spont=[];
doIsi=0;

if ~isempty(useAssign)
    spikes=filtspikes(spikes,0,'assigns',useAssign);
end

for i=1:length(useLEDcond)
    spikes=makeTempField(spikes,'led',useLEDcond);
    temp(i,:)=spikes.temp;
    temp1(i,:)=spikes.sweeps.temp;
    spikes.temp=[];
    spikes.sweeps.temp=[];
end
spikes.temp=sum(temp,1)>=1;
spikes.sweeps.temp=sum(temp1,1)>=1;
spikes=filtspikes(spikes,0,'temp',1,'stimcond',usestimcond);
% spikes=filtspikes(spikes,0,'led',useLEDcond,'stimcond',usestimcond);
s=spikes.spiketimes;

fre=[];
avSpec=[];
fre2=[];
avSpec2=[];
freSpont=[];
avSpecSpont=[];
if doCrossCorr==1
    if showPSTH==1
        figure(); 
    end
%     ax=axes; 
    ax=[];
    [h,ax,~,~,~,x,y]=psthForCycle_noShow(spikes,2,ax,0,stimWindow(2)-baseWindow(1));
    psthEv.x=x(x>=stimWindow(1) & x<=stimWindow(2))-x(find(x>=stimWindow(1),1,'first'));
    psthEv.y=y(x>=stimWindow(1) & x<=stimWindow(2));
    psthSp.x=x(x>=baseWindow(1) & x<baseWindow(2))-x(find(x>=baseWindow(1),1,'first'));
    psthSp.y=y(x>=baseWindow(1) & x<baseWindow(2));
    stim.x=psthEv.x;
    stim.y=sin(2*pi*floor(useLEDcond).*stim.x);
    [cc,cc_lags]=xcorr(stim.y-mean(stim.y),psthEv.y-mean(psthEv.y),[],'coeff'); % first input lags second input
    if length(stim.y(1:floor(end/2)))<length(psthSp.y)
        psthSp.y=psthSp.y(1:length(stim.y(1:floor(end/2))));
    end
    [ccSpont,cc_lagsSpont]=xcorr(stim.y(1:floor(end/2))-mean(stim.y(1:floor(end/2))),psthSp.y-mean(psthSp.y),[],'coeff'); 
    % neural response must lag stimulus
    % only consider as valid lags > 0
    new_ccX=cc_lags.*(stim.x(2)-stim.x(1));
    new_ccXSpont=cc_lagsSpont.*(stim.x(2)-stim.x(1));
%     disp(useLEDcond);
    mm=max(cc(new_ccX>0 & new_ccX<=1))-min(cc(new_ccX>0 & new_ccX<=1));
%     disp(mm);
    if showCrossCorr==1
        figure();
        plot(new_ccX,cc);
        title('Evoked');
        hold on; 
        figure(); 
        plot(new_ccXSpont,ccSpont);
        title('Spont');
    end
    sampRate=new_ccX(2)-new_ccX(1);
    spontSampRate=new_ccXSpont(2)-new_ccXSpont(1);
    if doFFT==1
        [fre,avSpec]=makePowerSpectrum(cc(new_ccX>0 & new_ccX<=1),1/sampRate);
        if showFigs==1
            figure(); 
            plot(fre,avSpec);
            title('Crosscorr Power Spec');
        end
        [freSpont,avSpecSpont]=makePowerSpectrum(ccSpont(new_ccXSpont>0 & new_ccXSpont<=1),1/spontSampRate);
    end
end

spiketimes=s(s>=stimWindow(1) & s<=stimWindow(2));
if length(spiketimes)>1
    if doIsi==1
        bins=round(1000*maxlag/corr_bin_size);
        isi=diff(spiketimes);
        isi=isi(isi<=maxlag);
        [n,x]=hist(isi*1000,linspace(0,1000*maxlag,bins));
        ymax=max(n)+1;
        lags=x;
        cross=n;
    else
%         maxlag=stimWindow(2)-stimWindow(1);
%         maxlag=0.1;
        maxlag=0.04;
        spiketimes=sort(spikes.unwrapped_times(s>=stimWindow(1) & s<=stimWindow(2)));        
        [cross,lags]=pxcorr(spiketimes,spiketimes,round(1000/corr_bin_size),maxlag);
        cross(lags==0)=0;
        ymax=max(cross)+1;
    end
    if doFFT==1
        
        params.Fs=1/(lags(2)-lags(1));
%         params.tapers=[2 3];
%         params.tapers=[1 1];
%         params.tapers=[5 maxlag*2 0];
        params.tapers=[9 maxlag*2 0];
        [avSpec2,fre2]=mtspectrumpb(cross,params);
        
%         d=lags(2)-lags(1);
%         [fre2,avSpec2]=makePowerSpectrum(cross,1/d);
        autofre=fre2(4:end);
        autoavSpec=avSpec2(4:end);
        forfigs.evoked_autocorr_power.x=autofre;
        forfigs.evoked_autocorr_power.y=autoavSpec;
        if showFigs==1
            figure(); 
            plot(autofre,autoavSpec);
            title('Autocorr Power Spec');
            xlim([30 80]);  
        end
    end
else
    cross=0;
    lags=0;
end

f1=figure();
hold on;
% for i=0:(1/floor(useLEDcond))*1000:maxlag*1000
%     if i==0
%     else
%         line([i i],[0 ymax],'Color','red');
%         line([-i -i],[0 ymax],'Color','red');
%     end
% end
% bb=bar(lags*1000,cross,1.0);
plot(lags*1000,cross);
forfigs.evoked_autocorr.x=lags*1000;
forfigs.evoked_autocorr.y=cross;
evoked.x=lags*1000;
evoked.y=cross;
% set(bb,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
% set(gca,'XLim',maxlag*1000*[-1 1]);
set(gca,'XLim',1000*[-0.1 0.1]);
set(gca,'YLim',[0 ymax]);
xlabel('Time lag (ms)');
if doIsi==1
    ylabel('Count');
else
    ylabel('Autocorrelation (Hz)');
end
title('Evoked');
if showFigs==0
    close(f1);
end

freSpont2=[];
avSpecSpont2=[];
spiketimes=s(s>=baseWindow(1) & s<=baseWindow(2));
if length(spiketimes)>1
    spiketimes=sort(spikes.unwrapped_times(s>=baseWindow(1) & s<=baseWindow(2)));        
    [cross,lags]=pxcorr(spiketimes,spiketimes,round(1000/corr_bin_size),maxlag);
    cross(lags==0)=0;
    ymax=max(cross)+1;
    if doFFT==1
        params.Fs=1/(lags(2)-lags(1));
        params.tapers=[0.9 maxlag*2 0];
        [avSpecSpont2,freSpont2]=mtspectrumpb(cross,params);
        
%         d=lags(2)-lags(1);
%         [freSpont2,avSpecSpont2]=makePowerSpectrum(cross,1/d);
        autofre=freSpont2(4:end);
        autoavSpec=avSpecSpont2(4:end);
        if showFigs==1
%             figure(); 
%             plot(autofre,autoavSpec);
%             title('Autocorr Power Spec - Spont');
        end
    end
else
    cross=0;
    lags=0;
end
% f2=figure();
% hold on;
% for i=0:(1/floor(useLEDcond))*1000:maxlag*1000
%     if i==0
%     else
%         line([i i],[0 ymax],'Color','red');
%         line([-i -i],[0 ymax],'Color','red');
%     end
% end
% bb=bar(lags*1000,cross,1.0);
% figure(); 
% plot(lags*1000,cross);
% title('Spont');
spont.x=lags*1000;
spont.y=cross;
% set(gca,'XLim',1000*[-0.1 0.1]);
% set(gca,'YLim',[0 ymax]);
% xlabel('Time lag (ms)');
% set(bb,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
% set(gca,'XLim',maxlag*1000*[-1 1]);
% set(gca,'YLim',[0 ymax]);
% xlabel('Time lag (ms)');
% ylabel('Autocorrelation (Hz)');
% title('Spontaneous');
% if showFigs==0
%     close(f2);
% end