function [outspikes,out]=matchAndMergeUnits(spikes1,spikes2,bestUnits1,bestUnits2,ITI,saveDir)

rpv_threshold=0.0015; % in ms
% sameUnitRPV_thresh=0.06;
sameUnitRPV_thresh=0.1;


if isempty(bestUnits1)
    bestUnits1=unique(spikes1.assigns);
end
if isempty(bestUnits2)
    bestUnits2=unique(spikes2.assigns);
end

f=figure();
movegui(f,'south');
disp(length(bestUnits1)*length(bestUnits2));
matchedPairs=[];
for i=1:length(bestUnits1)
    spikes_unit1=filtspikes(spikes1,0,'assigns',bestUnits1(i));
    unit1_wvfms=reshape(nanmean(spikes_unit1.waveforms,1),size(spikes_unit1.waveforms,2),4);
    for j=1:length(bestUnits2)
        % If overlap in sorted fileInd
        if any(ismember(spikes1.fileInd,spikes2.fileInd))
            overlapFileInd=unique(spikes1.fileInd(ismember(spikes1.fileInd,spikes2.fileInd)));
            % Check rpv in overlapping fileInd
            spikes_unit2=filtspikes(spikes2,0,'assigns',bestUnits2(j));
            sub1=filtspikes(spikes_unit1,0,'fileInd',overlapFileInd);
            sub2=filtspikes(spikes_unit2,0,'fileInd',overlapFileInd);
            sub1.trials=sub1.trials-min(sub1.trials)+1;
            sub2.trials=sub2.trials-min(sub2.trials)+1;
            sub1.unwrapped_times=(sub1.trials-1).*ITI+sub1.spiketimes;
            sub2.unwrapped_times=(sub2.trials-1).*ITI+sub2.spiketimes;
            if isempty(sub1.led) | isempty(sub2.led)
                continue
            end
            togSpikes=concatSpikes_shiftAssigns(sub1,sub2);
            togSpikes.unwrapped_times=sort(togSpikes.unwrapped_times);
            rpv=sum(diff(togSpikes.unwrapped_times)<rpv_threshold);
            out.rpvFraction=rpv./length(togSpikes.spiketimes);
            
            % Check whether # rpv greater than sameUnitRPV_thresh % of total spikes
            passes_rpv=rpv>=sameUnitRPV_thresh*length(togSpikes.spiketimes);
            if passes_rpv==0
                continue
            end
            disp(out.rpvFraction);
        end
        
        % Test whether waveforms match
        unit2_wvfms=reshape(nanmean(spikes_unit2.waveforms,1),size(spikes_unit2.waveforms,2),4);
        plot([unit1_wvfms(:,1); unit1_wvfms(:,2); unit1_wvfms(:,3); unit1_wvfms(:,4)],'Color','k');
        hold on;
        plot([unit2_wvfms(:,1); unit2_wvfms(:,2); unit2_wvfms(:,3); unit2_wvfms(:,4)],'Color','r');
        hold off;
        doesScale=questdlg('Do waveforms match?');
        if strcmp(doesScale,'No')
            passes_scale=0;
            continue
        elseif strcmp(doesScale,'Yes')
            passes_scale=1;
        elseif strcmp(doesScale,'Cancel')
            return
        end
        out.unit1wvfm=[unit1_wvfms(:,1); unit1_wvfms(:,2); unit1_wvfms(:,3); unit1_wvfms(:,4)];
        out.unit2wvfm=[unit2_wvfms(:,1); unit2_wvfms(:,2); unit2_wvfms(:,3); unit2_wvfms(:,4)];
        
        % Matched units
        matchedPairs=[matchedPairs; bestUnits1(i) bestUnits2(j)];
        disp('matched units');
    end
end

% Remove overlap
if any(ismember(spikes1.fileInd,spikes2.fileInd))
    f=unique(spikes2.fileInd);
    sub2=filtspikes(spikes2,0,'fileInd',f(~ismember(f,overlapFileInd)));
    spikes2=sub2;    
end

% Prepare trials
spikes1.trials=spikes1.trials-min(spikes1.trials)+1;
spikes1.sweeps.trials=spikes1.sweeps.trials-min(spikes1.sweeps.trials)+1;
spikes2.trials=spikes2.trials-min(spikes2.trials)+1;
spikes2.sweeps.trials=spikes2.sweeps.trials-min(spikes2.sweeps.trials)+1;
spikes2.trials=spikes2.trials+max(spikes1.trials);
spikes2.sweeps.trials=spikes2.sweeps.trials+max(spikes1.sweeps.trials);

% Make new conjoint spikes structure
spikes1.oldassigns=spikes1.assigns;
spikes1.assigns=spikes1.assigns+10000;
spikes2.oldassigns=spikes2.assigns;
spikes2.assigns=spikes2.assigns+10000;
for i=1:size(matchedPairs,1)
    spikes2.assigns(ismember(spikes2.oldassigns,matchedPairs(i,2)))=matchedPairs(i,1);    
end

outspikes=concatSpikes_noShiftAssigns(spikes1,spikes2);
% Rename assigns
[~,~,rnk]=unique(outspikes.assigns);
outspikes.assigns=rnk;

save([saveDir '\' 'outspikes.mat'],'outspikes');
save([saveDir '\' 'out.mat'],'out');
save([saveDir '\' 'matchedPairs.mat'],'matchedPairs');

end

function allTogether=concatSpikes_noShiftAssigns(spikes1,spikes2)

allTogether.led=[spikes1.led spikes2.led];
allTogether.stimcond=[spikes1.stimcond spikes2.stimcond];
try
    allTogether.waveforms=[spikes1.waveforms; spikes2.waveforms];
catch
    disp('hey');
end
allTogether.spiketimes=[spikes1.spiketimes spikes2.spiketimes];
allTogether.info.detect.event_channel=[spikes1.info.detect.event_channel; spikes2.info.detect.event_channel];
allTogether.trials=[spikes1.trials spikes2.trials];
allTogether.unwrapped_times=[spikes1.unwrapped_times spikes2.unwrapped_times];
allTogether.fileInd=[spikes1.fileInd spikes2.fileInd];
allTogether.trigger=[spikes1.trigger spikes2.trigger];
allTogether.time=[spikes1.time spikes2.time];
allTogether.assigns=[spikes1.assigns spikes2.assigns];

allTogether.sweeps.fileInd=[spikes1.sweeps.fileInd spikes2.sweeps.fileInd];
allTogether.sweeps.trials=[spikes1.sweeps.trials spikes2.sweeps.trials];
allTogether.sweeps.trigger=[spikes1.sweeps.trigger spikes2.sweeps.trigger];
allTogether.sweeps.stimcond=[spikes1.sweeps.stimcond spikes2.sweeps.stimcond];
allTogether.sweeps.led=[spikes1.sweeps.led spikes2.sweeps.led];

end
            
        