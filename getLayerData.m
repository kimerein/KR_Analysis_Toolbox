function getLayerData(spikes,recalc,fileInd)

dataDir='W:\Analysis Computer\Layers MUA\Awake\Mawake49\';
T='T4';

if recalc==1
    [~,evch]=recalculate_event_channel(spikes);
    spikes.event_channel=evch;
else
    spikes.event_channel=spikes.info.detect.event_channel';
end

for i=1:4
    [xpoints,ypoints1,ypoints2]=scriptForComparingMUA(filtspikes(spikes,0,'event_channel',i),fileInd,[],[]);
    data.xpoints=xpoints;
    data.ypoints1=ypoints1;
    data.ypoints2=ypoints2;
    save([dataDir T '_evch' num2str(i) '_led0vs5.mat'],'data');
%     save([dataDir T '_evch' num2str(i) '_led26vs48.mat'],'data');
end
close all;