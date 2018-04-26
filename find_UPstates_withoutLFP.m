function [UPstates,mua_xpoints,mua_ypoints]=find_UPstates_withoutLFP(spikes)

useTrials=1:length(unique(spikes.trials));
considerWindow=[0 3];

MUAthresh=12.56; % Get from 250 ms histogram

% If want to use MUA to define UP states
a=unique(spikes.trials);
UPstates_MUA=cell(length(useTrials));
totalSpikes=0;
spikesInUp=0;
for i=1:length(useTrials)
    subSpikes=filtspikes(spikes,0,'trials',a(useTrials(i)));
    [~,~,~,xpoints,ypoints]=psth_wStdev_valuesOnly(subSpikes,250,0);
    mua_xpoints=xpoints;
    mua_ypoints{i}=ypoints;
    UPstates_MUA{i}=returnUPs(xpoints,ypoints,MUAthresh);
    spiketimes=subSpikes.spiketimes;
    totalSpikes=totalSpikes+sum(spiketimes>considerWindow(1) & spiketimes<considerWindow(2));
    currUPs=UPstates_MUA{i};
    for j=1:size(currUPs,1)
       spikesInUp=spikesInUp+sum(spiketimes>currUPs(j,1) & spiketimes<currUPs(j,2) & spiketimes>considerWindow(1) & spiketimes<considerWindow(2));
    end
end
disp('totalSpikes');
disp(totalSpikes);
disp('spikesInUp');
disp(spikesInUp);
startTime=xpoints(1);
endTime=xpoints(end);

% Just use MUA UP states
UPstates=UPstates_MUA;