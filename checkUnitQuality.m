function [units_areBigEnough,useAssigns]=checkUnitQuality(spikes,useAssigns,ITI,detectionThreshold)

% rpv_threshold=0.0015; % in ms
rpv_threshold=0.001; % in s
% rpv_threshold=0; % in ms

if isempty(useAssigns)
    useAssigns=unique(spikes.assigns);
end

spikes.unwrapped_times=(spikes.trials-1).*ITI+spikes.spiketimes;

evCh=spikes.info.detect.event_channel;
f=figure();
movegui(f,'south');
units_areBigEnough=zeros(1,length(useAssigns));
for i=1:length(useAssigns)
    subspikes=filtspikes(spikes,0,'assigns',useAssigns(i));
    
    % Get # rpv
    rpv=sum(diff(subspikes.unwrapped_times)<rpv_threshold);
    rpvFraction=rpv./length(subspikes.spiketimes);
    % Check whether # rpv fewer than 1% of total spikes
    passes_rpv=rpv<=0.01*length(subspikes.spiketimes);
    if passes_rpv==0
        disp('failed rpv test');
        disp(useAssigns(i));
        disp(rpvFraction);
        continue
    end
    
    bestEvCh=mode(evCh(spikes.assigns==useAssigns(i)));
    wvfmAmps=reshape(nanmin(subspikes.waveforms(:,:,bestEvCh),[],2),size(subspikes.waveforms,1),1);
    [n,xout]=hist(wvfmAmps,500);
    plot(xout,n,'Color','k');
    hold on;
   if ~isempty(detectionThreshold)
       line([detectionThreshold detectionThreshold],[min(n) max(n)],'Color','r');
   elseif isfield(spikes.info.detect,'thresh')
       line([spikes.info.detect.thresh(bestEvCh) spikes.info.detect.thresh(bestEvCh)],[min(n) max(n)],'Color','r');
   end
   hold off;
   disp(useAssigns(i));
    isBigEnough=questdlg('Is unit big enough?');
    if strcmp(isBigEnough,'No')
        units_areBigEnough(i)=0;
    elseif strcmp(isBigEnough,'Yes')
        units_areBigEnough(i)=1;
    elseif strcmp(doesScale,'Cancel')
        return
    end
end
        
    