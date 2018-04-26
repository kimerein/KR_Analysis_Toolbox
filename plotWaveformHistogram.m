function allPtoT=plotWaveformHistogram(spikes,Fs,parameter)

if strcmp(parameter,'peakToTroughTime')
    allPtoT=zeros(size(spikes.waveforms,1),1);
    for i=1:size(spikes.waveforms,1)
        currEventCh=spikes.info.detect.event_channel(i);
        [peak,peakInd]=findpeaks(double(-spikes.waveforms(i,:,currEventCh)),'SORTSTR','descend','NPEAKS',1);
        if length(spikes.waveforms(i,:,currEventCh))-peakInd<2
            troughInd=length(spikes.waveforms(i,:,currEventCh));
        else
            [trough,troughInd]=findpeaks(double(spikes.waveforms(i,peakInd:end,currEventCh)),'SORTSTR','descend','NPEAKS',1);
        end
        if isempty(troughInd)
            troughInd=length(spikes.waveforms(i,:,currEventCh));
        else
            troughInd(1)=troughInd(1)+peakInd-1;
        end
        allPtoT(i)=troughInd(1)-peakInd(1);
    end
    allPtoT=allPtoT/Fs;
    
%     figure();
%     n=2;
%     plot(spikes.waveforms(n,:,spikes.info.detect.event_channel(n)));
%     currEventCh=spikes.info.detect.event_channel(n);
%     [peak,peakInd]=findpeaks(double(-spikes.waveforms(n,:,currEventCh)),'SORTSTR','descend','NPEAKS',1);
%     [trough,troughInd]=findpeaks(double(spikes.waveforms(n,peakInd:end,currEventCh)),'SORTSTR','descend','NPEAKS',1);
%     if isempty(troughInd)
%         troughInd=length(spikes.waveforms(n,:,currEventCh));
%     else
%         troughInd(1)=troughInd(1)+peakInd-1;
%     end
%     hold on;
%     plot([peakInd peakInd+0.02],[-peak -peak-0.02],'Color','red');
%     plot([troughInd troughInd+0.02],[trough trough-0.02],'Color','blue'); 
%     
%     figure();
%     plot(-spikes.waveforms(1,:,spikes.info.detect.event_channel(1)));
%     
    figure();
    hist(allPtoT,300);
    
    figure(); 
    c=histc(allPtoT,0:0.0001:0.001);
    bar(0:0.0001:0.001,c);
    
    figure();
    [n,xout]=hist(allPtoT,50);
    bar(xout,n);
end

if strcmp(parameter,'half-width')
    allPtoT=zeros(size(spikes.waveforms,1),1);
    for i=1:size(spikes.waveforms,1)
        currEventCh=spikes.info.detect.event_channel(i);
        shift1=double(-spikes.waveforms(i,:,currEventCh))-min(double(-spikes.waveforms(i,:,currEventCh)));
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
        allPtoT(i)=point2-point1;
        if i==1
            figure();
            plot(spikes.waveforms(i,:,spikes.info.detect.event_channel(i)));
            hold on;
            plot([point1 point1+0.02],[spikes.waveforms(i,point1,spikes.info.detect.event_channel(i)) spikes.waveforms(i,point1,spikes.info.detect.event_channel(i))-0.02],'Color','red');
            plot([point2 point2+0.02],[spikes.waveforms(i,point2,spikes.info.detect.event_channel(i)) spikes.waveforms(i,point2,spikes.info.detect.event_channel(i))-0.02],'Color','blue');
        end
    end
    allPtoT=allPtoT/Fs;

    figure();
    hist(allPtoT,300);
    
    figure(); 
    c=histc(allPtoT,0:0.0001:0.001);
    bar(0:0.0001:0.001,c);
    
    figure();
    [n,xout]=hist(allPtoT,50);
    bar(xout,n);
end