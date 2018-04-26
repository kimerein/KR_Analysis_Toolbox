function [thinSpikes,regSpikes,thinUnitsAssigns,regUnitsAssigns]=separateThinAndRegSpikingUnits(spikes,allPtoT)

if isempty(allPtoT)
    allPtoT=zeros(size(spikes.waveforms,1),1);
    for i=1:size(spikes.waveforms,1)
        currEventCh=spikes.info.detect.event_channel(i);
        [peak,peakInd]=findpeaks(double(-spikes.waveforms(i,:,currEventCh)),'SORTSTR','descend','NPEAKS',1);
        [trough,troughInd]=findpeaks(double(spikes.waveforms(i,peakInd:end,currEventCh)),'SORTSTR','descend','NPEAKS',1);
        if isempty(troughInd)
            troughInd=length(spikes.waveforms(i,:,currEventCh));
        else
            troughInd(1)=troughInd(1)+peakInd-1;
        end
        allPtoT(i)=troughInd(1)-peakInd(1);
    end
    allPtoT=allPtoT/Fs;
end
    
a=unique(spikes.assigns);
unitAvs=zeros(length(a),1);
for j=1:length(a)
    unitAvs(j)=mean(allPtoT(spikes.assigns==a(j)));
end

thinUnitsAssigns=[];
regUnitsAssigns=[];
for j=1:length(a)
    if unitAvs(j)<=7*10^-4
        thinUnitsAssigns=[thinUnitsAssigns; a(j)];
        disp(unitAvs(j));
        disp('thin');
    else
        regUnitsAssigns=[regUnitsAssigns; a(j)];
        disp(unitAvs(j));
        disp('reg');
    end
end

thinSpikes=filtspikes(spikes,0,'assigns',thinUnitsAssigns');
regSpikes=filtspikes(spikes,0,'assigns',regUnitsAssigns');

scriptForComparingMUA(thinSpikes);
scriptForComparingMUA(regSpikes);
