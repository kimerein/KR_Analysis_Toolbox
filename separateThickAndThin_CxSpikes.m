
thresh=3.2*10^-4;
ledConds=[3 6 7];

thinSpikes.spiketimes=[spikes.spiketimes(newHalfWidths<thresh) spikes2.spiketimes(oldHalfWidths<thresh) T5234_spikes_120206_files45to80.spiketimes(PVhalfWidths<thresh)];
thinSpikes.led2=[spikes.led(newHalfWidths<thresh) spikes2.led2(oldHalfWidths<thresh) T5234_spikes_120206_files45to80.led(PVhalfWidths<thresh)];

thickSpikes.spiketimes=[spikes.spiketimes(newHalfWidths>=thresh) spikes2.spiketimes(oldHalfWidths>=thresh) T5234_spikes_120206_files45to80.spiketimes(PVhalfWidths>=thresh)];
thickSpikes.led2=[spikes.led(newHalfWidths>=thresh) spikes2.led2(oldHalfWidths>=thresh) T5234_spikes_120206_files45to80.led(PVhalfWidths>=thresh)];
% thickSpikes.spiketimes=spikes.spiketimes(allHalfWidths>=thresh);
% thickSpikes.led2=spikes.led(allHalfWidths>=thresh);

[xpoints_spontSil,ypoints_spontSilThin]=psthWithOnlySpiketimes(thinSpikes.spiketimes(ismember(thinSpikes.led2,ledConds)),0.1,6,32000);
[xpoints_spontSil,ypoints_spontSilThick]=psthWithOnlySpiketimes(thickSpikes.spiketimes(ismember(thickSpikes.led2,ledConds)),0.1,6,32000);

ypoints_spontSilThin=ypoints_spontSilThin/length(ledConds);
ypoints_spontSilThick=ypoints_spontSilThick/length(ledConds);

%
ledConds=[1 2 4];

[xpoints_spontSil,ypoints_spontNoThin]=psthWithOnlySpiketimes(thinSpikes.spiketimes(ismember(thinSpikes.led2,ledConds)),0.1,6,32000);
[xpoints_spontSil,ypoints_spontNoThick]=psthWithOnlySpiketimes(thickSpikes.spiketimes(ismember(thickSpikes.led2,ledConds)),0.1,6,32000);

ypoints_spontNoThin=ypoints_spontNoThin/length(ledConds);
ypoints_spontNoThick=ypoints_spontNoThick/length(ledConds);