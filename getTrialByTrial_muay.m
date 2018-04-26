function [mua_y,xpoints]=getTrialByTrial_muay(spikes,fileInds)

spikes=filtspikes(spikes,0,'fileInd',fileInds);

trials=unique(spikes.trials);
mua_y=cell(1,length(trials));
for i=1:length(trials)
    disp(i);
    [xpoints,ypoints1]=scriptForComparingMUA(spikes,[],trials(i),trials(i));
    mua_y{i}=ypoints1;
end




