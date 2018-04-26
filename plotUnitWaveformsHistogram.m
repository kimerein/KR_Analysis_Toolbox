function allPtoT=plotUnitWaveformsHistogram(spikes,Fs,parameter)

if strcmp(parameter,'peakToTroughTime')
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
    
    a=unique(spikes.assigns);
    unitAvs=zeros(length(a),1);
    for j=1:length(a)
        unitAvs(j)=mean(allPtoT(spikes.assigns==a(j)));
    end
    
    
    figure();
    [n,xout]=hist(unitAvs,50);
    bar(xout,n);
end