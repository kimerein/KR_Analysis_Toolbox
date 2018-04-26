function [units_areBigEnough,useAssigns]=plotWaveformAmpDistribution(spikes,useAssigns)


if isempty(useAssigns)
    useAssigns=unique(spikes.assigns);
end

evCh=spikes.info.detect.event_channel;
f=figure();
movegui(f,'south');
units_areBigEnough=nan(1,length(useAssigns));
for i=1:length(useAssigns)
    subspikes=filtspikes(spikes,0,'assigns',useAssigns(i));
    bestEvCh=mode(evCh(spikes.assigns==useAssigns(i)));
    wvfmAmps=reshape(nanmin(subspikes.waveforms(:,:,bestEvCh),[],2),size(subspikes.waveforms,1),1);
    [n,xout]=hist(wvfmAmps,500);
    plot(xout,n,'Color','k');
    hold on;
   if isfield(spikes.info.detect,'thresh')
       line([spikes.info.detect.thresh(bestEvCh) spikes.info.detect.thresh(bestEvCh)],[min(n) max(n)],'Color','r');
   end
   hold off;
    isBigEnough=questdlg('Is unit big enough?');
    if strcmp(isBigEnough,'No')
        units_areBigEnough(i)=0;
    elseif strcmp(isBigEnough,'Yes')
        units_areBigEnough(i)=1;
    elseif strcmp(doesScale,'Cancel')
        return
    end
end
        
    