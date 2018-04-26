function [s,r]=calculateUnitsSuppressionByLED(spikes)

s=[];
a=unique(spikes.assigns);
for i=1:length(a)
    noLED=filtspikes(spikes,0,'assigns',a(i),'led',0);
    LED=filtspikes(spikes,0,'assigns',a(i),'led',5);
    noLEDspikes=sum(((noLED.spiketimes>1) + (noLED.spiketimes<2))-1);
    LEDspikes=sum(((LED.spiketimes>1) + (LED.spiketimes<2))-1);
    % Have to normalize by number of trials and length of interval
    noLEDspikes=noLEDspikes/length(unique(noLED.trials));
    LEDspikes=LEDspikes/length(unique(LED.trials));
    s=[s; noLEDspikes-LEDspikes];
end
    
r=[];
for i=1:length(a)
    noLED=filtspikes(spikes,0,'assigns',a(i),'led',0);
    LED=filtspikes(spikes,0,'assigns',a(i),'led',5);
    noLEDspikes=sum(((noLED.spiketimes>2) + (noLED.spiketimes<2.5))-1);
    LEDspikes=sum(((LED.spiketimes>2) + (LED.spiketimes<2.5))-1);
    % Have to normalize by number of trials and length of interval
    noLEDspikes=noLEDspikes/(length(unique(noLED.trials))/0.5);
    LEDspikes=LEDspikes/(length(unique(LED.trials))/0.5);
    r=[r; LEDspikes-noLEDspikes];
end

    