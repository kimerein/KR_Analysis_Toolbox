function [evoked,spont,useLEDcond,maxlag,fre,avSpec]=plot_unit_autocorr_earlyEpoch(spikes,useAssign,useLEDcond)

showFigs=0;
doCrossCorr=1;
showCrossCorr=0;
doFFT=1;
showPSTH=0;
corr_bin_size=2; % in ms
maxlag=2; % in s
% useLEDcond=60;
usestimcond=1:9;
stimWindow=[1 3];
earlyEpoch=[1 1.5];
baseWindow=[0 1];

if ~isempty(useAssign)
    spikes=filtspikes(spikes,0,'assigns',useAssign);
end

spikes=filtspikes(spikes,0,'led',useLEDcond,'stimcond',usestimcond);
s=spikes.spiketimes;

fre=[];
avSpec=[];
if doCrossCorr==1
    if showPSTH==1
        figure(); 
    end
%     ax=axes; 
    ax=[];
    [h,ax,~,~,~,x,y]=psthForCycle_noShow(spikes,2,ax,0,stimWindow(2)-baseWindow(1));
    psthEv.x=x(x>=earlyEpoch(1) & x<=earlyEpoch(2))-x(find(x>=earlyEpoch(1),1,'first'));
    psthEv.y=y(x>=earlyEpoch(1) & x<=earlyEpoch(2));
    psthSp.x=x(x>=baseWindow(1) & x<baseWindow(2))-x(find(x>=baseWindow(1),1,'first'));
    psthSp.y=y(x>=baseWindow(1) & x<baseWindow(2));
    stim.x=psthEv.x;
    stim.y=sin(2*pi*useLEDcond.*stim.x);
    [cc,cc_lags]=xcorr(stim.y-mean(stim.y),psthEv.y-mean(psthEv.y),[],'coeff'); % first input lags second input
    % neural response must lag stimulus
    % only consider as valid lags > 0
    new_ccX=cc_lags.*(stim.x(2)-stim.x(1));
%     disp(useLEDcond);
    mm=max(cc(new_ccX>0 & new_ccX<=1))-min(cc(new_ccX>0 & new_ccX<=1));
%     disp(mm);
    if showCrossCorr==1
        figure();
        plot(new_ccX,cc);
    end
    sampRate=new_ccX(2)-new_ccX(1);
    if doFFT==1
        [fre,avSpec]=makePowerSpectrum(cc(new_ccX>0 & new_ccX<=1),1/sampRate);
    end
end

spiketimes=s(s>=stimWindow(1) & s<=stimWindow(2));
if length(spiketimes)>1
    [cross,lags]=pxcorr(spiketimes,spiketimes,round(1000/corr_bin_size),maxlag);
%     if doFFT==1
%         d=lags(2)-lags(1);
%         [fre2,avSpec2]=makePowerSpectrum(cross,1/d);
%         figure(); 
%         plot(fre2,avSpec2);
%         title('Autocorr Power Spec');
%     end
else
    cross=0;
    lags=0;
end
cross(find(lags==0))=0;
ymax=max(cross)+1;
f1=figure();
hold on;
for i=0:(1/useLEDcond)*1000:maxlag*1000
    if i==0
    else
        line([i i],[0 ymax],'Color','red');
        line([-i -i],[0 ymax],'Color','red');
    end
end
bb=bar(lags*1000,cross,1.0);
evoked.x=lags*1000;
evoked.y=cross;
set(bb,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
set(gca,'XLim',maxlag*1000*[-1 1]);
set(gca,'YLim',[0 ymax]);
xlabel('Time lag (ms)');
ylabel('Autocorrelation (Hz)');
title('Evoked');
if showFigs==0
    close(f1);
end

spiketimes=s(s>=baseWindow(1) & s<=baseWindow(2));
if length(spiketimes)>1
    [cross,lags]=pxcorr(spiketimes,spiketimes,round(1000/corr_bin_size),maxlag);
else
    cross=0;
    lags=0;
end
cross(find(lags==0))=0;
ymax=max(cross)+1;
f2=figure();
hold on;
for i=0:(1/useLEDcond)*1000:maxlag*1000
    if i==0
    else
        line([i i],[0 ymax],'Color','red');
        line([-i -i],[0 ymax],'Color','red');
    end
end
bb=bar(lags*1000,cross,1.0);
spont.x=lags*1000;
spont.y=cross;
set(bb,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
set(gca,'XLim',maxlag*1000*[-1 1]);
set(gca,'YLim',[0 ymax]);
xlabel('Time lag (ms)');
ylabel('Autocorrelation (Hz)');
title('Spontaneous');
if showFigs==0
    close(f2);
end