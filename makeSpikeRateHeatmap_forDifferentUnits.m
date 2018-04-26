function [xpoints,spikeRates,order]=makeSpikeRateHeatmap_forDifferentUnits(spikes,ledConds,order)

spikes=filtspikes(spikes,0,'led',ledConds);

a=unique(spikes.assigns);
for i=1:length(a)
    disp(i);
    [n,centers,edges,xpoints,ypoints]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(i)),50,1);
    spikeRates(i,:)=ypoints;
end

if isempty(order)
    [sortedFRs,order]=sort(mean(spikeRates,2),'descend');
    spikeRates=spikeRates(order,:);
else
    spikeRates=spikeRates(order,:);
end

figure();
% heatmap(spikeRates);
leg=min(spikeRates(1:end)):(max(spikeRates(1:end))-min(spikeRates(1:end)))/size(spikeRates,2):max(spikeRates(1:end));
imagesc([spikeRates; leg(1:length(spikeRates(1,:)))]);
