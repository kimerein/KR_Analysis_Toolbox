function [thinSpikes,regSpikes,thinUnitsAssigns,regUnitsAssigns,unitAvs]=separateThinAndRegSpikingUnits_byAvWaveform(spikes,allPtoT)

% if isempty(allPtoT)
%     allPtoT=zeros(size(spikes.waveforms,1),1);
%     for i=1:size(spikes.waveforms,1)
%         currEventCh=spikes.info.detect.event_channel(i);
%         [peak,peakInd]=findpeaks(double(-spikes.waveforms(i,:,currEventCh)),'SORTSTR','descend','NPEAKS',1);
%         [trough,troughInd]=findpeaks(double(spikes.waveforms(i,peakInd:end,currEventCh)),'SORTSTR','descend','NPEAKS',1);
%         if isempty(troughInd)
%             troughInd=length(spikes.waveforms(i,:,currEventCh));
%         else
%             troughInd(1)=troughInd(1)+peakInd-1;
%         end
%         allPtoT(i)=troughInd(1)-peakInd(1);
%     end
%     allPtoT=allPtoT/Fs;
% end
    
a=unique(spikes.assigns);

avWaveforms=zeros(length(a),size(spikes.waveforms,2),size(spikes.waveforms,3));
evChsForAssigns={};
for j=1:length(a)
    avWaveforms(j,:,:)=mean(spikes.waveforms(spikes.assigns==a(j),:,:),1);
    evChsForAssigns{j}=spikes.info.detect.event_channel(spikes.assigns==a(j));
end
    
unitAvs=zeros(length(a),1);
for j=1:length(a)
    currEventCh=mode(evChsForAssigns{j});
    shift1=double(-avWaveforms(j,:,currEventCh))-min(double(-avWaveforms(j,:,currEventCh)));
    [peak,peakInd]=findpeaks(shift1,'SORTSTR','descend','NPEAKS',1);
    halfAmp=(peak-shift1(1))/2;
    shift2=shift1-halfAmp;
    point1=find(shift2(peakInd:-1:1)<0,1,'first');
    if isempty(point1)
        point1=1;
    end
    point1=peakInd-point1+1;
    point2=find(shift2(peakInd:end)<0,1,'first');
    if isempty(point2)
        point2=length(shift2(peakInd:end));
    end
    point2=peakInd+point2-1;
    unitAvs(j)=point2-point1;
    %unitAvs(j)=halfAmp*2;
end
unitAvs=unitAvs/32000;

thinUnitsAssigns=[];
regUnitsAssigns=[];
for j=1:length(a)
    if unitAvs(j)<=3*10^-4
%     if unitAvs(j)<=0.2
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

scriptForComparingMUA(thinSpikes,[],[],[]);
scriptForComparingMUA(regSpikes,[],[],[]);
