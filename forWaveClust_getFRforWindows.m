function [sup,reb,vis,spont]=forWaveClust_getFRforWindows(spikes,ledVal,noLedVal)

sup=[];
a=unique(spikes.assigns);
w(1)=1;
w(2)=1.4;
% for i=1:length(a)
%     noLED=filtspikes(spikes,0,'assigns',a(i),'led',noLedVal);
%     LED=filtspikes(spikes,0,'assigns',a(i),'led',ledVal);
%     noLEDspikes=sum( (noLED.spiketimes>w(1)) + (noLED.spiketimes<w(2)) );
%     LEDspikes=sum((LED.spiketimes>w(1)) + (LED.spiketimes<w(2)));
% %     noLEDspikes=sum(((noLED.spiketimes>w(1)) + (noLED.spiketimes<w(2)))-1);
% %     LEDspikes=sum(((LED.spiketimes>w(1)) + (LED.spiketimes<w(2)))-1);
%     % Have to normalize by number of trials and length of interval
%     noLEDspikes=noLEDspikes/(length(unique(noLED.trials))*(w(2)-w(1)));
%     LEDspikes=LEDspikes/(length(unique(LED.trials))*(w(2)-w(1)));
%     sup=[sup; noLEDspikes-LEDspikes];
% end
for i=1:length(a)
    noLED=filtspikes(spikes,0,'assigns',a(i),'led',noLedVal);
    [m,s,n]=calcMeanAndStdDuringWindow(noLED,[w(1) w(2)]);
    [m2,s2,n2]=calcMeanAndStdDuringWindow(noLED,[0 1]);
    LED=filtspikes(spikes,0,'assigns',a(i),'led',ledVal);
    [m,s,LED_n]=calcMeanAndStdDuringWindow(LED,[w(1) w(2)]);
    [m2,s2,LED_n2]=calcMeanAndStdDuringWindow(LED,[0 1]);
    sup=[sup; mean(n-n2)-mean(LED_n-LED_n2)];
end
    
% reb=[];
% w(1)=4.3;
% w(2)=5;
% for i=1:length(a)
%     noLED=filtspikes(spikes,0,'assigns',a(i),'led',noLedVal);
%     LED=filtspikes(spikes,0,'assigns',a(i),'led',ledVal);
%     noLEDspikes=sum(((noLED.spiketimes>w(1)) + (noLED.spiketimes<w(2)))-1);
%     LEDspikes=sum(((LED.spiketimes>w(1)) + (LED.spiketimes<w(2)))-1);
%     % Have to normalize by number of trials and length of interval
%     noLEDspikes=noLEDspikes/(length(unique(noLED.trials))*(w(2)-w(1)));
%     LEDspikes=LEDspikes/(length(unique(LED.trials))*(w(2)-w(1)));
%     reb=[reb; noLEDspikes-LEDspikes];
% end

reb=[];
w(1)=0.54;
w(2)=1;
for i=1:length(a)
    noLED=filtspikes(spikes,0,'assigns',a(i),'led',noLedVal);
    LED=filtspikes(spikes,0,'assigns',a(i),'led',ledVal);
    noLEDspikes=sum((noLED.spiketimes>w(1)) + (noLED.spiketimes<w(2)));
    LEDspikes=sum((LED.spiketimes>w(1)) + (LED.spiketimes<w(2)));
    % Have to normalize by number of trials and length of interval
    noLEDspikes=noLEDspikes/(length(unique(noLED.trials))*(w(2)-w(1)));
    LEDspikes=LEDspikes/(length(unique(LED.trials))*(w(2)-w(1)));
%     reb=[reb; LEDspikes-noLEDspikes];
    reb=[reb; noLEDspikes-LEDspikes];
end

vis=[];
w(1)=1;
w(2)=1.4;
% for i=1:length(a)
%     noLED=filtspikes(spikes,0,'assigns',a(i),'led',noLedVal);
%     noLEDspikes=sum((noLED.spiketimes>w(1)) + (noLED.spiketimes<w(2)));
%     % Have to normalize by number of trials and length of interval
%     noLEDspikes=noLEDspikes/(length(unique(noLED.trials))*(w(2)-w(1)));
%     vis=[vis; noLEDspikes];
% end    
for i=1:length(a)
    noLED=filtspikes(spikes,0,'assigns',a(i),'led',noLedVal);
    [m,s,n]=calcMeanAndStdDuringWindow(noLED,[w(1) w(2)]);
    [m2,s2,n2]=calcMeanAndStdDuringWindow(noLED,[0 1]);
    vis=[vis; mean(n-n2)];
%     LED=filtspikes(spikes,0,'assigns',a(i),'led',ledVal);
%     [m,s,LED_n]=calcMeanAndStdDuringWindow(LED,[w(1) w(2)]);
%     [m2,s2,LED_n2]=calcMeanAndStdDuringWindow(LED,[0 1]);
%     vis=[vis; mean(n-n2)-mean(LED_n-LED_n2)];
end

spont=[];
w(1)=0;
w(2)=0.46;
for i=1:length(a)
    noLED=filtspikes(spikes,0,'assigns',a(i),'led',noLedVal);
    LED=filtspikes(spikes,0,'assigns',a(i),'led',ledVal);
    noLEDspikes=sum((noLED.spiketimes>w(1)) + (noLED.spiketimes<w(2)));
    LEDspikes=sum((LED.spiketimes>w(1)) + (LED.spiketimes<w(2)));
    % Have to normalize by number of trials and length of interval
    noLEDspikes=noLEDspikes/(length(unique(noLED.trials))*(w(2)-w(1)));
    LEDspikes=LEDspikes/(length(unique(LED.trials))*(w(2)-w(1))); 
    spont=[spont; mean([LEDspikes noLEDspikes])];
end