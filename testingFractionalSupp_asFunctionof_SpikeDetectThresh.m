function testingFractionalSupp_asFunctionof_SpikeDetectThresh(spikes)

spikeDetectThreshs=[0.01 0.02 0.03 0.04 0.05];

% Cut down spikes using spikeDetectThresh
for i=1:length(spikeDetectThreshs)
    thresh=spikeDetectThreshs(i);
    acrossChThresh=spikes.info.detect.thresh+thresh;
    amps=max(spikes.waveforms,1)-min(spikes.waveforms,1);
    spikes.keep=amps>
    spikes = makeTempField(spikes,'led',cond.values{n});
    cspikes(m,n) = filtspikes(spikes,0,'stimcond',stim.code{m},'temp',1);
    
    dtSteps=0.001:0.001:0.043;
    spontAc=[1.05*10^4 8620 7250 6000 5100 4200 3400 2750 2250 1800 1500 1200 950 750 650 500 380 330 225 185 125 100 90 70 55 50 50 40 35 20 18 18 12 10 10 7 5 2.5 2.5 2.5 2.5 2.5 0];
    
    spontSil_LEDconds=[2 7];
    spontNo_LEDconds=[1 3];
    dtSteps=0.001:0.001:0.045;
    minWvfrms=min(wvfrms.a,[],2);
    supp_spontAc=zeros(length(dtSteps),1);
    noSupp_spontAc=zeros(length(dtSteps),1);
    fracSupp_spontAc=zeros(length(dtSteps),1);
    for i=1:length(dtSteps)
        th=infoDOTdetect.a.thresh(5)-dtSteps(i);
        subSpikes=spikes1.spiketimes(minWvfrms(:,:,5)<th & infoDOTdetect.a.event_channel==5 & ismember(spikes1.led,spontSil_LEDconds));
        [xpoints,ypoints]=psthWithOnlySpiketimes(subSpikes,0.01,6,32000);
        close all;
        supp_spontAc(i)=mean(ypoints(xpoints>1.7 & xpoints<1.9));
        subSpikes=spikes1.spiketimes(minWvfrms(:,:,5)<th & infoDOTdetect.a.event_channel==5 & ismember(spikes1.led,spontSil_LEDconds));
        [xpoints,ypoints]=psthWithOnlySpiketimes(subSpikes,0.01,6,32000);
        close all;
        noSupp_spontAc(i)=mean(ypoints(xpoints>1.7 & xpoints<1.9));
        fracSupp_spontAc(i)=supp_spontAc(i)/noSupp_spontAc(i);
    end
    figure();
    plot(dtSteps,fracSupp_spontAc);
        
        
        