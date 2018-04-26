detectThresh=3.2*10^-4;
spontSil_LEDconds=[2 7];
spontNo_LEDconds=[1 3];
dtSteps=-0.04:0.001:0.06;
minWvfrms=min(wvfrms.a,[],2);
supp_spontAc=zeros(length(dtSteps),1);
noSupp_spontAc=zeros(length(dtSteps),1);
fracSupp_spontAc=zeros(length(dtSteps),1);
for i=1:length(dtSteps)
%     th=infoDOTdetect.a.thresh-dtSteps(i);
%     numT=length(unique(spikes1.trials((newHalfWidths>detectThresh)' & min(minWvfrms(:,:,:),[],3)<th & infoDOTdetect.a.event_channel==5 & ismember(spikes1.led,spontSil_LEDconds)')));
%     subSpikes=spikes1.spiketimes((newHalfWidths>detectThresh)' & min(minWvfrms(:,:,:),[],3)<th & infoDOTdetect.a.event_channel==5 & ismember(spikes1.led,spontSil_LEDconds)');
    detectCutoffs=zeros(size(min(minWvfrms(:,:,:),[],3)));
    for j=1:16
        th=infoDOTdetect.a.thresh(j)-dtSteps(i);
        detectCutoffs(min(minWvfrms(:,:,j),[],3)<th & infoDOTdetect.a.event_channel==j)=1;
    end
    numT=length(unique(spikes1.trials((newHalfWidths>=detectThresh)' & detectCutoffs & ismember(spikes1.led,spontSil_LEDconds)')));
    subSpikes=spikes1.spiketimes((newHalfWidths>=detectThresh)' & detectCutoffs & ismember(spikes1.led,spontSil_LEDconds)');
    [xpoints,ypoints]=psthWithOnlySpiketimes(subSpikes,0.01,6,32000);
    if i~=1 close all; end
    supp_spontAc(i)=mean(ypoints(xpoints>1.7 & xpoints<1.9))/(numT*0.2);
%     numT=length(unique(spikes1.trials((newHalfWidths>detectThresh)' & min(minWvfrms(:,:,:),[],3)<th & infoDOTdetect.a.event_channel==5 & ismember(spikes1.led,spontNo_LEDconds)')));
%     subSpikes=spikes1.spiketimes((newHalfWidths>detectThresh)' & min(minWvfrms(:,:,:),[],3)<th & infoDOTdetect.a.event_channel==5 & ismember(spikes1.led,spontNo_LEDconds)');
    numT=length(unique(spikes1.trials((newHalfWidths>=detectThresh)' & detectCutoffs & ismember(spikes1.led,spontNo_LEDconds)')));
    subSpikes=spikes1.spiketimes((newHalfWidths>=detectThresh)' & detectCutoffs & ismember(spikes1.led,spontNo_LEDconds)');
    [xpoints,ypoints]=psthWithOnlySpiketimes(subSpikes,0.01,6,32000);
    if i~=1 close all; end
    noSupp_spontAc(i)=mean(ypoints(xpoints>1.7 & xpoints<1.9))/(numT*0.2);
%     fracSupp_spontAc(i)=noSupp_spontAc(i)./supp_spontAc(i);
    fracSupp_spontAc(i)=supp_spontAc(i)./noSupp_spontAc(i);
%     fracSupp_spontAc(i)=(noSupp_spontAc(i)-supp_spontAc(i))./noSupp_spontAc(i);
end
figure();
plot(dtSteps,fracSupp_spontAc);
figure();
plot(dtSteps,noSupp_spontAc,'Color','k');
hold on;
plot(dtSteps,supp_spontAc,'Color','r');