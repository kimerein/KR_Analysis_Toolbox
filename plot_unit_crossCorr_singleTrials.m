function [evoked,spont,useLEDcond,maxlag,allFre,allAvSpec,fre2,avSpec2,allFreSpont,allAvSpecSpont,freSpont2,avSpecSpont2]=plot_unit_crossCorr_singleTrials(spikes,useAssign,useLEDcond,useFreq)

showFigs=0;
doCrossCorr=1;
showCrossCorr=0;
doFFT=1;
showPSTH=0;
corr_bin_size=2; % in ms
maxlag=2; % in s
% useLEDcond=60;
usestimcond=1:128;
stimWindow=[1 3];
% stimWindow=[1.3 3];
% stimWindow=[1.3 3];
baseWindow=[0 1];
evoked=[];
spont=[];
allFre=[];
allAvSpec=[];
fre2=[];
avSpec2=[];
allFreSpont=[];
allAvSpecSpont=[];
freSpont2=[];
avSpecSpont2=[];
doPower=1;

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
    t=unique(spikes.trials);
    if doPower==1
        for i=1:length(t)
            [h,ax,~,~,~,x,y]=psthForCycle_noShow(filtspikes(spikes,0,'trials',t),2,ax,0,stimWindow(2)-baseWindow(1));
            psthEv.x=x(x>=stimWindow(1) & x<=stimWindow(2))-x(find(x>=stimWindow(1),1,'first'));
            psthEv.y=y(x>=stimWindow(1) & x<=stimWindow(2));
            psthEv.y=psthEv.y-mean(psthEv.y);
            psthSp.x=x(x>=baseWindow(1) & x<baseWindow(2))-x(find(x>=baseWindow(1),1,'first'));
            psthSp.y=y(x>=baseWindow(1) & x<baseWindow(2));
            if i==1
                allEvTrials=zeros(length(t),length(psthEv.y));
                allEvTrials(i,:)=psthEv.y;
            else
                allEvTrials(i,:)=psthEv.y;
            end
        end
        sampRate=psthEv.x(2)-psthEv.x(1);
        params.Fs=1/sampRate;
        params.tapers=[0.9 stimWindow(2)-stimWindow(1) 0];
        [S1,f1]=mtspectrumpb(allEvTrials',params);
        allFre=f1;
        allAvSpec=S1;
    else
        for i=1:length(t)
            [h,ax,~,~,~,x,y]=psthForCycle_noShow(filtspikes(spikes,0,'trials',t),2,ax,0,stimWindow(2)-baseWindow(1));
            psthEv.x=x(x>=stimWindow(1) & x<=stimWindow(2))-x(find(x>=stimWindow(1),1,'first'));
            psthEv.y=y(x>=stimWindow(1) & x<=stimWindow(2));
            stim.x=psthEv.x;
            stim.y=sin(2*pi*floor(useFreq).*stim.x);
            psthSp.x=x(x>=baseWindow(1) & x<baseWindow(2))-x(find(x>=baseWindow(1),1,'first'));
            psthSp.y=y(x>=baseWindow(1) & x<baseWindow(2));
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
            if i==1
                allFre=zeros(length(t),length(fre));
                allAvSpec=zeros(length(t),length(avSpec));
                allFreSpont=zeros(length(t),length(freSpont));
                allAvSpecSpont=zeros(length(t),length(avSpecSpont));
                allFre(i,:)=fre;
                allAvSpec(i,:)=avSpec;
                allFreSpont(i,:)=freSpont;
                allAvSpecSpont(i,:)=avSpecSpont;
            end
        end
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