function [evoked,spont,useLEDcond,maxlag,fre,avSpec,fre2,avSpec2,freSpont,avSpecSpont,freSpont2,avSpecSpont2]=plot_LFP_crossCorr(LFPbySweep,useLEDcond,useFreq,ledBySweep,stimBySweep)

showFigs=0;
doCrossCorr=1;
showCrossCorr=0;
doFFT=1;
showPSTH=0;
corr_bin_size=2; % in ms
maxlag=2; % in s
% useLEDcond=60;
usestimcond=1:128;
% stimWindow=[1 3];
stimWindow=[1.3 3];
baseWindow=[0 1];
evoked=[];
spont=[];
trialDuration=4;

lfp=LFPbySweep(ismember(ledBySweep,useLEDcond) & ismember(stimBySweep,usestimcond),:);

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
    x=linspace(0,trialDuration,size(lfp,2));
    y=mean(lfp,1);
    psthEv.x=x(x>=stimWindow(1) & x<=stimWindow(2))-x(find(x>=stimWindow(1),1,'first'));
    psthEv.y=y(x>=stimWindow(1) & x<=stimWindow(2));
    psthSp.x=x(x>=baseWindow(1) & x<baseWindow(2))-x(find(x>=baseWindow(1),1,'first'));
    psthSp.y=y(x>=baseWindow(1) & x<baseWindow(2));
    stim.x=psthEv.x;
    stim.y=sin(2*pi*floor(useFreq).*stim.x);
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

% spiketimes=s(s>=stimWindow(1) & s<=stimWindow(2));
% if length(spiketimes)>1
%     [cross,lags]=pxcorr(spiketimes,spiketimes,round(1000/corr_bin_size),maxlag);
%     if doFFT==1
%         d=lags(2)-lags(1);
%         [fre2,avSpec2]=makePowerSpectrum(cross,1/d);
%         autofre=fre2(4:end);
%         autoavSpec=avSpec2(4:end);
%         if showFigs==1
%             figure(); 
%             plot(autofre,autoavSpec);
%             title('Autocorr Power Spec');
%         end
%     end
% else
%     cross=0;
%     lags=0;
% end
% cross(find(lags==0))=0;
% ymax=max(cross)+1;
% if showFigs==1
%     f1=figure();
%     hold on;
%     for i=0:(1/floor(useFreq))*1000:maxlag*1000
%         if i==0
%         else
%             line([i i],[0 ymax],'Color','red');
%             line([-i -i],[0 ymax],'Color','red');
%         end
%     end
%     bb=bar(lags*1000,cross,1.0);
% end
% evoked.x=lags*1000;
% evoked.y=cross;
% if showFigs==1
%     set(bb,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
%     set(gca,'XLim',maxlag*1000*[-1 1]);
%     set(gca,'YLim',[0 ymax]);
%     xlabel('Time lag (ms)');
%     ylabel('Autocorrelation (Hz)');
%     title('Evoked');
% end
% % if showFigs==0
% %     close(f1);
% % end

% freSpont2=[];
% avSpecSpont2=[];
% spiketimes=s(s>=baseWindow(1) & s<=baseWindow(2));
% if length(spiketimes)>1
%     [cross,lags]=pxcorr(spiketimes,spiketimes,round(1000/corr_bin_size),maxlag);
%     if doFFT==1
%         d=lags(2)-lags(1);
%         [freSpont2,avSpecSpont2]=makePowerSpectrum(cross,1/d);
%         autofre=freSpont2(4:end);
%         autoavSpec=avSpecSpont2(4:end);
%         if showFigs==1
%             figure(); 
%             plot(autofre,autoavSpec);
%             title('Autocorr Power Spec - Spont');
%         end
%     end
% else
%     cross=0;
%     lags=0;
% end
% cross(find(lags==0))=0;
% ymax=max(cross)+1;
% if showFigs==1
%     f2=figure();
%     hold on;
%     for i=0:(1/floor(useFreq))*1000:maxlag*1000
%         if i==0
%         else
%             line([i i],[0 ymax],'Color','red');
%             line([-i -i],[0 ymax],'Color','red');
%         end
%     end
% end
% spont.x=lags*1000;
% spont.y=cross;
% if showFigs==1
%     bb=bar(lags*1000,cross,1.0);
%     set(bb,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
%     set(gca,'XLim',maxlag*1000*[-1 1]);
%     set(gca,'YLim',[0 ymax]);
%     xlabel('Time lag (ms)');
%     ylabel('Autocorrelation (Hz)');
%     title('Spontaneous');
% end
% if showFigs==0
%     close(f2);
% end