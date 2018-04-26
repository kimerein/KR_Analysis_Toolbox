function [UPstates,powerRatio,mua_xpoints,mua_ypoints]=find_UPstates_with_LFP(LFPbySweep,LFP_Fs,lowerBand,upperBand,spikes)

global dontShowSpecgram
dontShowSpecgram=1;

dontGetMUAUPs=1;
showExample=1;
exampleTrial=1;
UPthresh=9.2;
% UPthresh=4;
useTrials=1:size(LFPbySweep,1);
% useTrials=unique(spikes.trials);
useTrials=useTrials-useTrials(1)+1;
considerWindow=[0 5];

% MUAthresh=30.3; % Get from 250 ms histogram
MUAthresh=49.6; % Get from 250 ms histogram

powerRatio=[];

if showExample
    dontShowSpecgram=0;
    getPowerIntegrals(LFPbySweep,exampleTrial,LFP_Fs,lowerBand,upperBand,1,UPthresh,spikes,MUAthresh);
end
UPstates_LFP=cell(length(useTrials),1);
dontShowSpecgram=1;
powerRatio=cell(length(useTrials),1);
for i=1:length(useTrials)
    [UPstates_LFP{i},powerRatio{i}]=getPowerIntegrals(LFPbySweep,useTrials(i),LFP_Fs,lowerBand,upperBand,0,UPthresh,[],MUAthresh);
end

% If also want to use MUA to define UP states
a=unique(spikes.trials);
% if length(a)~=size(LFPbySweep,1)
%     disp('spikes # of trials and size of LFPbySweep do not match!');
%     return
% end
UPstates_MUA=cell(1,length(useTrials));
totalSpikes=0;
spikesInUp=0;
for i=1:length(useTrials)
    %     subSpikes=filtspikes(spikes,0,'trials',a(useTrials(i)));
    subSpikes=filtspikes(spikes,0,'trials',useTrials(i));
    if isempty(subSpikes)
        continue
    end
    %     subSpikes=filtspikes(spikes,0,'trials',useTrials(i));
    %     [~,~,~,xpoints,ypoints]=psth_wStdev_valuesOnly(subSpikes,250,0);
    [~,~,~,xpoints,ypoints]=psth_wStdev_valuesOnly(subSpikes,250,0);
    mua_xpoints=xpoints;
    mua_ypoints{i}=ypoints;
    if dontGetMUAUPs==0
        UPstates_MUA{i}=returnUPs(xpoints,ypoints,MUAthresh);
        spiketimes=subSpikes.spiketimes;
        totalSpikes=totalSpikes+sum(spiketimes>considerWindow(1) & spiketimes<considerWindow(2));
        currUPs=UPstates_MUA{i};
        for j=1:size(currUPs,1)
            spikesInUp=spikesInUp+sum(spiketimes>currUPs(j,1) & spiketimes<currUPs(j,2) & spiketimes>considerWindow(1) & spiketimes<considerWindow(2));
        end
    end
end
if dontGetMUAUPs==0
    disp('totalSpikes');
    disp(totalSpikes);
    disp('spikesInUp');
    disp(spikesInUp);
end
startTime=xpoints(1);
endTime=xpoints(end);


% Reconcile UP states identified by LFP power ratio and MUA
% UPstates=cell(length(useTrials));
% for i=1:length(useTrials)
%     times=startTime:0.001:endTime;
%     timesInUP=zeros(1,length(times));
%     currTrialUPs=UPstates_LFP{i};
%     for j=1:size(currTrialUPs,1)
% %         timesInUP(times>=currTrialUPs(j,1) & times<=currTrialUPs(j,2))=1;
%         timesInUP(times>=currTrialUPs(j,1) & times<=currTrialUPs(j,2))=0.4;
%     end
%     currTrialUPs=UPstates_MUA{i};
%     for j=1:size(currTrialUPs,1)
% %         timesInUP(times>=currTrialUPs(j,1) & times<=currTrialUPs(j,2))=1;
%         timesInUP(times>=currTrialUPs(j,1) & times<=currTrialUPs(j,2))=timesInUP(times>=currTrialUPs(j,1) & times<=currTrialUPs(j,2))+0.4;
%     end
%     UPstates{i}=returnUPs(times,timesInUP,0.5);
% end
% UPstates=cell(length(useTrials));

% Just use MUA UP states
% UPstates=UPstates_MUA;
UPstates=UPstates_LFP;