function [unit_wvfms,unit_halfWidths,unit_depths,useAssigns]=getUnitWaveformAndDepth(spikes)

useAssigns=unique(spikes.assigns);
evCh=spikes.info.detect.event_channel;

if length(spikes.info.detect.event_channel)~=length(spikes.spiketimes)
    % Concatenated spikes struct; thus, evCh is wrong
    % Need to recalculate
    maxCh=nan(1,length(useAssigns));
    for i=1:length(useAssigns)
        ch1_wvfm=reshape(nanmean(spikes.waveforms(spikes.assigns==useAssigns(i),:,1),1),1,size(spikes.waveforms,2));
        ch2_wvfm=reshape(nanmean(spikes.waveforms(spikes.assigns==useAssigns(i),:,2),1),1,size(spikes.waveforms,2));
        ch3_wvfm=reshape(nanmean(spikes.waveforms(spikes.assigns==useAssigns(i),:,3),1),1,size(spikes.waveforms,2));
        ch4_wvfm=reshape(nanmean(spikes.waveforms(spikes.assigns==useAssigns(i),:,4),1),1,size(spikes.waveforms,2));
        ch1_min=min(ch1_wvfm);
        ch2_min=min(ch2_wvfm);
        ch3_min=min(ch3_wvfm);
        ch4_min=min(ch4_wvfm);
        [~,maxCh(i)]=min([ch1_min ch2_min ch3_min ch4_min]);
    end
    evCh=nan(size(spikes.assigns));
    for i=1:length(useAssigns)
        evCh(spikes.assigns==useAssigns(i))=maxCh(i);
    end    
    
    spikes.params.Fs=25000;
end

unit_wvfms=nan(length(useAssigns),size(spikes.waveforms,2));
unit_halfWidths=nan(length(useAssigns),1);

for i=1:length(useAssigns)
    bestEvCh=mode(evCh(spikes.assigns==useAssigns(i)));
    wvfms=reshape(spikes.waveforms(spikes.assigns==useAssigns(i),:,bestEvCh),sum(spikes.assigns==useAssigns(i)),size(spikes.waveforms,2));
    av_wvfm=nanmean(wvfms,1);

    halfWidth=classifyUnitWaveformWidth(av_wvfm,0.22*10^-3,spikes.params.Fs);
    
    unit_wvfms(i,:)=av_wvfm;
    unit_halfWidths(i)=halfWidth;
end

unit_depths=nan(length(useAssigns),1);
evchs=orderUnitsbyDepth_getEvChs(spikes,useAssigns);

for i=1:length(useAssigns)
    curr=evchs(i,:);
    [m,ind]=max(curr);
    newcurr=curr;
    newcurr(ind)=-10;
    [m2,ind2]=max(newcurr);
    newcurr2=newcurr;
    newcurr2(ind2)=-10;
    [m3,ind3]=max(newcurr2);
    if curr(ind)==curr(ind2)
        unit_depths(i)=mean([ind ind2]);
    elseif curr(ind2)==curr(ind3)
        unit_depths(i)=ind;
    else
        a=curr(ind);
        b=curr(ind2);
        y=(100*b)/(a+b);
        x=100-y;
        if ind>ind2
            unit_depths(i)=ind-(y/100);
        else
            unit_depths(i)=ind+(y/100);
        end
    end
end
    
    
    
    
    
    
    
    
    
    
   
    
    
    