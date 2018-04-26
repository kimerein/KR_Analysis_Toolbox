function plotUnitRawPSTH(spikes)

freq1=2;
freq2=20;

a=unique(spikes.assigns);

backSpikes=spikes;
for i=1:length(a)
    currA=a(i);
    disp(i);
    spikes=filtspikes(backSpikes,0,'assigns',currA);
    spikesLED1=sortSpikesForLED(spikes,freq1+0.05);
    spikesLED2=sortSpikesForLED(spikes,freq2+0.05);
    ax=[];
    [~,~,~,~,~,x_freq1,y_freq1(i,:)]=psthForCycle_noShow(spikesLED1,10,ax,0,4);
    [~,~,~,~,~,x_freq2,y_freq2(i,:)]=psthForCycle_noShow(spikesLED2,10,ax,0,4);
end

figure();
runningOffset=0;
curr1=y_freq1;
curr2=y_freq2;
for i=1:length(a)
    curr1(i,:)=curr1(i,:)-min(curr1(i,:))+runningOffset;
    plot(x_freq1,curr1(i,:),'Color','k');
    hold on;
    curr2(i,:)=curr2(i,:)-min(curr2(i,:))+runningOffset;
    plot(x_freq2+max(x_freq1)+0.5,curr2(i,:),'Color','b');
    runningOffset=max([curr1(i,:) curr2(i,:)]);
end
title('units PSTH for 2 Hz (black) and 20 Hz (blue)');

end

function spikes=sortSpikesForLED(spikes,useLEDcond)

for i=1:length(useLEDcond)
    spikes=makeTempField(spikes,'led',useLEDcond);
    temp(i,:)=spikes.temp;
    temp1(i,:)=spikes.sweeps.temp;
    spikes.temp=[];
    spikes.sweeps.temp=[];
end
spikes.temp=sum(temp,1)>=1;
spikes.sweeps.temp=sum(temp1,1)>=1;
spikes=filtspikes(spikes,0,'temp',1);
end