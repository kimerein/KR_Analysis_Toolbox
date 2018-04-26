function newSpikes = subsampleSpikes(spikes, oneInX)
%
% Randomly subsamples a spikes structure
% 
% spikes        spikes struct
% oneInX        will randomly choose one in X (oneInX) spikes from spikes
%               to include in newSpikes
%
% Created: 11/3/11 - KR

newSpikes=spikes;
fields=fieldnames(spikes);
nSpikes=length(spikes.spiketimes);
selectThisMany=round(nSpikes/oneInX);
allSpikes=1:nSpikes;
selected=allSpikes(randi(length(allSpikes),[1 selectThisMany])); 
for i=1:length(fields)
    if isnumeric(spikes.(fields{i}))
        newSpikes.(fields{i})=spikes.(fields{i})(selected);
    end
end
        