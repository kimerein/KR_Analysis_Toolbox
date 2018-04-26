function [evoked,spont,useLEDcond,maxlag,fre,avSpec,fre2,avSpec2,freSpont,avSpecSpont,freSpont2,avSpecSpont2]=plot_unit_crossCorr(spikes,useAssign,useLEDcond,useFreq)

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
% stimWindow=[0.5 1];
% stimWindow=[1.4 3];
baseWindow=[0 1];
evoked=[];
spont=[];
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
% spikes=filtspikes(spikes,0,'temp',1,'stimcond',usestimcond);
spikes=filtspikes(spikes,0,'temp',1);
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
    [h,ax,~,~,~,x,y]=psthForCycle_noShow(spikes,2,ax,0,4);
%     [h,ax,~,~,~,x,y]=psthForCycle_noShow(spikes,2,ax,0,stimWindow(2)-baseWindow(1));
    psthEv.x=x(x>=stimWindow(1) & x<=stimWindow(2))-x(find(x>=stimWindow(1),1,'first'));
    if isempty(y)
        y=zeros(size(x));
        disp('Not Enough Trials');
    end
    psthEv.y=y(x>=stimWindow(1) & x<=stimWindow(2));
    % Uncomment next 3 for cross-corr
%     psthEv.y=psthEv.y-min(psthEv.y);
%     psthEv.y=psthEv.y./max(psthEv.y);
%     psthEv.y=psthEv.y-mean(psthEv.y)/2;
    psthSp.x=x(x>=baseWindow(1) & x<baseWindow(2))-x(find(x>=baseWindow(1),1,'first'));
    psthSp.y=y(x>=baseWindow(1) & x<baseWindow(2));
    % Use next line for raw amp. spectrum
%     psthEv.y=psthEv.y-mean(psthSp.y);
    psthSp.y=psthSp.y-mean(psthSp.y);
    stim.x=psthEv.x;
    stim.y=sin(2*pi*floor(useFreq).*stim.x);
%     [cc,cc_lags]=xcorr(stim.y-mean(stim.y),psthEv.y-mean(psthEv.y),[],'coeff'); % first input lags second input
    [cc,cc_lags]=xcorr(stim.y-mean(stim.y),psthEv.y-mean(psthEv.y),[],'none'); % first input lags second input
    if length(stim.y(1:floor(end/2)))<length(psthSp.y)
        psthSp.y=psthSp.y(1:length(stim.y(1:floor(end/2))));
    end
    [ccSpont,cc_lagsSpont]=xcorr(stim.y(1:floor(end/2))-mean(stim.y(1:floor(end/2))),psthSp.y-mean(psthSp.y),[],'none'); 
%     [ccSpont,cc_lagsSpont]=xcorr(stim.y(1:floor(end/2))-mean(stim.y(1:floor(end/2))),psthSp.y-mean(psthSp.y),[],'coeff'); 
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
        if doPower==0
            % Raw Amp Spectrum
            L=length(psthEv.x);
            Fs=1/0.002;
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            Y = fft(psthEv.y,NFFT)/L;
            fff = Fs/2*linspace(0,1,NFFT/2+1);
            fre=fff;
            avSpec=2*abs(Y(1:NFFT/2+1));
            useI2=find(fre>=useFreq,1,'first');
            useI=find(fre<=useFreq,1,'last');
            if abs(useI2-useFreq)<abs(useI-useFreq)
                avSpec(1:useI2-1)=0;
                avSpec(useI2+1:end)=0;
            else
                avSpec(1:useI-1)=0;
                avSpec(useI+1:end)=0;
            end
            %         avSpec(useI:useI2)=mean([avSpec(useI) avSpec(useI2)]);
            %         avSpec(1:useI-1)=0;
            %         avSpec(useI2+1:end)=0;
            %         if useFreq==1 & useLEDcond==1.05
            %             disp('hi');
            %         end
        else
            params.Fs=1/(psthEv.x(2)-psthEv.x(1));
            params.tapers=[0.9 stimWindow(2)-stimWindow(1) 0];
            [S1,f1]=mtspectrumpb(psthEv.y-mean(psthEv.y),params);
%             [S1,f1]=mtspectrumpb(cc(new_ccX>0 & new_ccX<=1),params); % RIGHT WAY
%             [S1,f1]=mtspectrumpb(cc(new_ccX>=0),params); % OK RIGHT WAY
            fre=f1;
            avSpec=S1;
%             [fre,avSpec]=makePowerSpectrum(cc(new_ccX>0 & new_ccX<=1),1/sampRate);
        end
        if showFigs==1
            figure(); 
            plot(fre,avSpec);
            title('Crosscorr Power Spec');
        end
        if doPower==0
            L=length(psthSp.x);
            Fs=1/0.002;
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            Y = fft(psthSp.y,NFFT)/L;
            fff = Fs/2*linspace(0,1,NFFT/2+1);
            freSpont=fff;
            avSpecSpont=2*abs(Y(1:NFFT/2+1));
            useI2=find(freSpont>=useFreq,1,'first');
            useI=find(freSpont<=useFreq,1,'last');
            if abs(useI2-useFreq)<abs(useI-useFreq)
                avSpecSpont(1:useI2-1)=0;
                avSpecSpont(useI2+1:end)=0;
            else
                avSpecSpont(1:useI-1)=0;
                avSpecSpont(useI+1:end)=0;
            end
            %         avSpecSpont(useI:useI2)=mean([avSpecSpont(useI) avSpecSpont(useI2)]);
            %         avSpecSpont(1:useI-1)=0;
            %         avSpecSpont(useI2+1:end)=0;
        else
            params.Fs=1/spontSampRate;
            params.tapers=[5 9];
            [S1,f1]=mtspectrumpb(ccSpont(new_ccXSpont>0 & new_ccXSpont<=1),params);
            freSpont=f1;
            avSpecSpont=S1;
%             [freSpont,avSpecSpont]=makePowerSpectrum(ccSpont(new_ccXSpont>0 & new_ccXSpont<=1),1/spontSampRate);
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