function [burstPairs]=checkForBursting(spikes,useAssigns,ITI,saveDir,bestUnits)

rpv_threshold=0.0015; % in ms
followfraction_threshold=0.05; % fraction of spikes that must follow
max_follow_interval=0.01; % in ms

% if ~isempty(burstPairsIn)
%     burstingUnits=combineBurstPairs(burstPairsIn);
%     save([saveDir '\burstingUnits.mat'],'burstingUnits');
%     return
% end

spikes.unwrapped_times=(spikes.trials-1).*ITI+spikes.spiketimes;

if isempty(useAssigns)
    useAssigns=unique(spikes.assigns);
end

if isempty(bestUnits)
    
f=figure();
movegui(f,'south');
burstPairs=[];
disp(length(useAssigns)*length(useAssigns));
for i=1:length(useAssigns)
        spikes_unit1=filtspikes(spikes,0,'assigns',useAssigns(i));
        unit1_wvfms=reshape(nanmean(spikes_unit1.waveforms,1),size(spikes_unit1.waveforms,2),4);
    for j=i+1:length(useAssigns)
        togSpikes=filtspikes(spikes,0,'assigns',[useAssigns(i) useAssigns(j)]);
        rpv=sum(diff(togSpikes.unwrapped_times)<rpv_threshold);
        out.rpvFraction=rpv./length(togSpikes.spiketimes);
        
        % Check whether # rpv fewer than 1% of total spikes
        passes_rpv=rpv<=0.01*length(togSpikes.spiketimes);
        if passes_rpv==0
            continue
        end
        
        % Test whether waveforms scale
        spikes_unit2=filtspikes(spikes,0,'assigns',useAssigns(j));
        unit2_wvfms=reshape(nanmean(spikes_unit2.waveforms,1),size(spikes_unit2.waveforms,2),4);
        plot([unit1_wvfms(:,1); unit1_wvfms(:,2); unit1_wvfms(:,3); unit1_wvfms(:,4)],'Color','k');
        hold on;
        plot([unit2_wvfms(:,1); unit2_wvfms(:,2); unit2_wvfms(:,3); unit2_wvfms(:,4)],'Color','r');
        hold off;
        doesScale=questdlg('Do waveforms scale?');
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
        
        % Test whether smaller waveform has fewer spikes
        unit1_min=min([unit1_wvfms(:,1); unit1_wvfms(:,2); unit1_wvfms(:,3); unit1_wvfms(:,4)]);
        unit2_min=min([unit2_wvfms(:,1); unit2_wvfms(:,2); unit2_wvfms(:,3); unit2_wvfms(:,4)]);
        if abs(unit1_min)<abs(unit2_min)  
            if length(spikes_unit1.spiketimes)<length(spikes_unit2.spiketimes)
                passes_fewerspikes=1;
            else
                passes_fewerspikes=0;
            end
        elseif abs(unit1_min)>abs(unit2_min)
            if length(spikes_unit1.spiketimes)>length(spikes_unit2.spiketimes)
                passes_fewerspikes=1;
            else
                passes_fewerspikes=0;
            end
        else
            passes_fewerspikes=0;
        end
        if passes_fewerspikes==0
            continue
        end
        disp('passed fewer spikes');
        
        % Test whether spikes in unit with smaller waveforms follow spikes
        % in larger unit
        togSpikes.isUnit=zeros(1,length(togSpikes.spiketimes));
        togSpikes.isSmallUnit=zeros(1,length(togSpikes.spiketimes));
        if abs(unit1_min)<abs(unit2_min)
            togSpikes.isUnit(togSpikes.assigns==useAssigns(i))=1;
            togSpikes.isUnit(togSpikes.assigns==useAssigns(j))=2;
            togSpikes.isSmallUnit(togSpikes.assigns==useAssigns(i))=1;
            togSpikes.isSmallUnit(togSpikes.assigns==useAssigns(j))=nan;
            numberSmallSpikes=length(spikes_unit1.spiketimes);
        else
            togSpikes.isUnit(togSpikes.assigns==useAssigns(j))=1;
            togSpikes.isUnit(togSpikes.assigns==useAssigns(i))=2;
            togSpikes.isSmallUnit(togSpikes.assigns==useAssigns(j))=1;
            togSpikes.isSmallUnit(togSpikes.assigns==useAssigns(i))=nan;
            numberSmallSpikes=length(spikes_unit2.spiketimes);
        end
        isi=diff(togSpikes.unwrapped_times);
        isBetweenUnits=diff(togSpikes.isUnit)==1 | ~isnan(diff(togSpikes.isSmallUnit));
        numberFollowSpikes=sum(isi<=max_follow_interval & logical(isBetweenUnits));
        fractionFollowSpikes=numberFollowSpikes/numberSmallSpikes;
        if fractionFollowSpikes>followfraction_threshold
            burstPairs=[burstPairs; useAssigns(i) useAssigns(j)];
        else
            continue
        end  
        disp('passed follow spikes');
        disp([useAssigns(i) useAssigns(j)]);
        out.isi_of_smallUnitSpikes=isi(isBetweenUnits==1);
        save([saveDir '\pair' num2str(useAssigns(i)) 'and' num2str(useAssigns(j)) '.mat'],'out');
    end
end

else
    
f=figure();
movegui(f,'south');
burstPairs=[];
disp(length(bestUnits)*length(useAssigns));
for i=1:length(bestUnits)
        spikes_unit1=filtspikes(spikes,0,'assigns',bestUnits(i));
        unit1_wvfms=reshape(nanmean(spikes_unit1.waveforms,1),size(spikes_unit1.waveforms,2),4);
    for j=1:length(useAssigns)
        togSpikes=filtspikes(spikes,0,'assigns',[bestUnits(i) useAssigns(j)]);
        rpv=sum(diff(togSpikes.unwrapped_times)<rpv_threshold);
        out.rpvFraction=rpv./length(togSpikes.spiketimes);
        
        % Check whether # rpv fewer than 1% of total spikes
        passes_rpv=rpv<=0.01*length(togSpikes.spiketimes);
        if passes_rpv==0
            continue
        end
        
        % Test whether waveforms scale
        spikes_unit2=filtspikes(spikes,0,'assigns',useAssigns(j));
        unit2_wvfms=reshape(nanmean(spikes_unit2.waveforms,1),size(spikes_unit2.waveforms,2),4);
        plot([unit1_wvfms(:,1); unit1_wvfms(:,2); unit1_wvfms(:,3); unit1_wvfms(:,4)],'Color','k');
        hold on;
        plot([unit2_wvfms(:,1); unit2_wvfms(:,2); unit2_wvfms(:,3); unit2_wvfms(:,4)],'Color','r');
        hold off;
        doesScale=questdlg('Do waveforms scale?');
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
        
        % Test whether smaller waveform has fewer spikes
        unit1_min=min([unit1_wvfms(:,1); unit1_wvfms(:,2); unit1_wvfms(:,3); unit1_wvfms(:,4)]);
        unit2_min=min([unit2_wvfms(:,1); unit2_wvfms(:,2); unit2_wvfms(:,3); unit2_wvfms(:,4)]);
        if abs(unit1_min)<abs(unit2_min)  
            if length(spikes_unit1.spiketimes)<length(spikes_unit2.spiketimes)
                passes_fewerspikes=1;
            else
                passes_fewerspikes=0;
            end
        elseif abs(unit1_min)>abs(unit2_min)
            if length(spikes_unit1.spiketimes)>length(spikes_unit2.spiketimes)
                passes_fewerspikes=1;
            else
                passes_fewerspikes=0;
            end
        else
            passes_fewerspikes=0;
        end
        if passes_fewerspikes==0
            continue
        end
        disp('passed fewer spikes');
        
        % Test whether spikes in unit with smaller waveforms follow spikes
        % in larger unit
        togSpikes.isUnit=zeros(1,length(togSpikes.spiketimes));
        togSpikes.isSmallUnit=zeros(1,length(togSpikes.spiketimes));
        if abs(unit1_min)<abs(unit2_min)
            togSpikes.isUnit(togSpikes.assigns==bestUnits(i))=1;
            togSpikes.isUnit(togSpikes.assigns==useAssigns(j))=2;
            togSpikes.isSmallUnit(togSpikes.assigns==bestUnits(i))=1;
            togSpikes.isSmallUnit(togSpikes.assigns==useAssigns(j))=nan;
            numberSmallSpikes=length(spikes_unit1.spiketimes);
        else
            togSpikes.isUnit(togSpikes.assigns==useAssigns(j))=1;
            togSpikes.isUnit(togSpikes.assigns==bestUnits(i))=2;
            togSpikes.isSmallUnit(togSpikes.assigns==useAssigns(j))=1;
            togSpikes.isSmallUnit(togSpikes.assigns==bestUnits(i))=nan;
            numberSmallSpikes=length(spikes_unit2.spiketimes);
        end
        isi=diff(togSpikes.unwrapped_times);
        isBetweenUnits=diff(togSpikes.isUnit)==1 | ~isnan(diff(togSpikes.isSmallUnit));
        numberFollowSpikes=sum(isi<=max_follow_interval & logical(isBetweenUnits));
        fractionFollowSpikes=numberFollowSpikes/numberSmallSpikes;
        if fractionFollowSpikes>followfraction_threshold
            burstPairs=[burstPairs; bestUnits(i) useAssigns(j)];
        else
            continue
        end  
        disp('passed follow spikes');
        disp([bestUnits(i) useAssigns(j)]);
        out.isi_of_smallUnitSpikes=isi(isBetweenUnits==1);
        save([saveDir '\pair' num2str(bestUnits(i)) 'and' num2str(useAssigns(j)) '.mat'],'out');
    end
end
    
end

save([saveDir '\burstPairs.mat'],'burstPairs');

% Combine all pairs
% burstingUnits=combineBurstPairs(burstPairs);
%     
% save([saveDir '\burstingUnits.mat'],'burstingUnits');


end

% function burstingUnits=combineBurstPairs(burstPairs)
% 
% % Combine all pairs
% allBurstUnits=burstPairs(1:end);
% burstingUnits=[];
% k=1;
% for i=1:length(allBurstUnits)
%     currUnit=allBurstUnits(i);
%     if isnan(currUnit)
%         continue
%     end
%     stillSearching=1;
%     pairedUnitSet=[];
%     saveInd=1;
%     while stillSearching==1 & saveInd<10000
%         takeRows=sum(ismember(burstPairs,currUnit),2);
%         pairsFromRows=burstPairs(takeRows>0,:);
%         uniqueBurstUnits=unique(pairsFromRows(1:end));
%         pairedUnitSet=[pairedUnitSet uniqueBurstUnits];
%     end
%     takeRows=sum(burstPairs==currUnit,2);
%     pairsFromRows=burstPairs(logical(takeRows),:);
%     uniqueBurstUnits=unique(pairsFromRows(1:end));
%     allBurstUnits(ismember(allBurstUnits,uniqueBurstUnits))=nan;
%     burstingUnits{k}=uniqueBurstUnits;
%     k=k+1;
% end
% 
% end

function allTogether=concatSpikes_sub(spikes1,spikes2)

allTogether.led=[spikes1.led spikes2.led];
allTogether.stimcond=[spikes1.stimcond spikes2.stimcond];
try
    allTogether.waveforms=[spikes1.waveforms; spikes2.waveforms];
catch
    disp('Waveform concatenation did not work.');
end
allTogether.spiketimes=[spikes1.spiketimes spikes2.spiketimes];
% Use the following line if DON'T want to shift event channels
% allTogether.info.detect.event_channel=[spikes1.info.detect.event_channel; spikes2.info.detect.event_channel];
% Else shift event channels
allTogether.info.detect.event_channel=[spikes1.info.detect.event_channel; spikes2.info.detect.event_channel+max(spikes1.info.detect.event_channel)];
allTogether.trials=[spikes1.trials spikes2.trials];
allTogether.unwrapped_times=[spikes1.unwrapped_times spikes2.unwrapped_times];
allTogether.fileInd=[spikes1.fileInd spikes2.fileInd];
allTogether.trigger=[spikes1.trigger spikes2.trigger];
allTogether.time=[spikes1.time spikes2.time];

allTogether.sweeps.fileInd=[spikes1.sweeps.fileInd spikes2.sweeps.fileInd];
allTogether.sweeps.trials=[spikes1.sweeps.trials spikes2.sweeps.trials];
allTogether.sweeps.trigger=[spikes1.sweeps.trigger spikes2.sweeps.trigger];
allTogether.sweeps.stimcond=[spikes1.sweeps.stimcond spikes2.sweeps.stimcond];
allTogether.sweeps.led=[spikes1.sweeps.led spikes2.sweeps.led];

[allTogether.trials,si]=sort(allTogether.trials);
allTogether.spiketimes=allTogether.spiketimes(si);
allTogether.led=allTogether.led(si);
allTogether.stimcond=allTogether.stimcond(si);
allTogether.waveforms=allTogether.waveforms(si,:);
allTogether.unwrapped_times=allTogether.unwrapped_times(si);
allTogether.fileInd=allTogether.fileInd(si);
allTogether.trigger=allTogether.trigger(si);
allTogether.time=allTogether.time(si);
allTogether.info.detect.event_channel=allTogether.info.detect.event_channel(si);
end