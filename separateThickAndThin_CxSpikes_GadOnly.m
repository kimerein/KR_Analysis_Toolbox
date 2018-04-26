function [rethinSpikes rethickSpikes xpoints_spontSil ypoints_spontNoThin ypoints_spontNoThick ypoints_spontSilThin ypoints_spontSilThick]=separateThickAndThin_CxSpikes_GadOnly(newHalfWidths,oldHalfWidths,PVhalfWidths,T5234_spikes_120206_files45to80,spikes1,spikes2)

thresh=3.2*10^-4;

Gad_sil=[6];
PV_sil=[3 8 9 10];
Gad_no=[1 2];
PV_no=[2 4];

a=spikes1.spiketimes(newHalfWidths<thresh & ismember(spikes1.led,Gad_sil)); % Gad 
b=spikes2.spiketimes(oldHalfWidths<thresh & ismember(spikes2.led2,Gad_sil)); % Gad
thinSpikes.spiketimes=[a b];

a=spikes1.trials(newHalfWidths<thresh & ismember(spikes1.led,Gad_sil)); % Gad
b=spikes2.trials(oldHalfWidths<thresh & ismember(spikes2.led2,Gad_sil)); % Gad 
thinSpikes.trials=[a b];

thickSpikes.spiketimes=[spikes1.spiketimes(newHalfWidths>=thresh & ismember(spikes1.led,Gad_sil))...
                    spikes2.spiketimes(oldHalfWidths>=thresh & ismember(spikes2.led2,Gad_sil))];

thickSpikes.trials=[spikes1.trials(newHalfWidths>=thresh & ismember(spikes1.led,Gad_sil))...
                    spikes2.trials(oldHalfWidths>=thresh & ismember(spikes2.led2,Gad_sil))];

[xpoints_spontSil,ypoints_spontSilThin]=psthWithOnlySpiketimes(thinSpikes.spiketimes,0.05,6,32000);
[xpoints_spontSil,ypoints_spontSilThick]=psthWithOnlySpiketimes(thickSpikes.spiketimes,0.05,6,32000);

% ypoints_spontSilThin=ypoints_spontSilThin/length([Gad_sil PV_sil]);
% ypoints_spontSilThick=ypoints_spontSilThick/length([Gad_sil PV_sil]);

rethinSpikes=thinSpikes;
rethickSpikes=thickSpikes;

%%%%%%%%%%%%%%%%%%%%%%%%%

thinSpikes.spiketimes=[spikes1.spiketimes(newHalfWidths<thresh & ismember(spikes1.led,Gad_no))...
                       spikes2.spiketimes(oldHalfWidths<thresh & ismember(spikes2.led2,Gad_no))];
thickSpikes.spiketimes=[spikes1.spiketimes(newHalfWidths>=thresh & ismember(spikes1.led,Gad_no))...
                        spikes2.spiketimes(oldHalfWidths>=thresh & ismember(spikes2.led2,Gad_no))];
thinSpikes.trials=[spikes1.trials(newHalfWidths<thresh & ismember(spikes1.led,Gad_no))...
                       spikes2.trials(oldHalfWidths<thresh & ismember(spikes2.led2,Gad_no))];
thickSpikes.trials=[spikes1.trials(newHalfWidths>=thresh & ismember(spikes1.led,Gad_no))...
                       spikes2.trials(oldHalfWidths>=thresh & ismember(spikes2.led2,Gad_no))];

% [xpoints_spontSil,ypoints_spontNoThin]=psthWithOnlySpiketimes(thinSpikes.spiketimes(ismember(thinSpikes.led2,ledConds)),0.1,6,32000);
% [xpoints_spontSil,ypoints_spontNoThick]=psthWithOnlySpiketimes(thickSpikes.spiketimes(ismember(thickSpikes.led2,ledConds)),0.1,6,32000);
[xpoints_spontSil,ypoints_spontNoThin]=psthWithOnlySpiketimes(thinSpikes.spiketimes,0.05,6,32000);
[xpoints_spontSil,ypoints_spontNoThick]=psthWithOnlySpiketimes(thickSpikes.spiketimes,0.05,6,32000);

% ypoints_spontNoThin=ypoints_spontNoThin/length([Gad_no PV_no]);
% ypoints_spontNoThick=ypoints_spontNoThick/length([Gad_no PV_no]);

% rethinSpikes=thinSpikes;
% rethickSpikes=thickSpikes;
