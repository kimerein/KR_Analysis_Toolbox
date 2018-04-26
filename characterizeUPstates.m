function [FRhist,durhist,numUPhist,FRs,durations,numUPs,trialAvFR,trialAvduration,trialdurationHist]=characterizeUPstates(UPstates,mua_xpoints,mua_ypoints)

correctUPend=0;

useUPsInTimeWindow=[0 3];
endsBy=1;
UPdurationThresh=0.2;

% UPstates=UPstates(:,1);

% Get UP state distributions for 
% Firing Rate
% Duration
% Number of UP states per trial
FRs=[];
durations=[];
numUPs=zeros(1,length(UPstates));
trialAvFR=zeros(1,length(UPstates));
trialAvduration=zeros(1,length(UPstates));
trialdurationHist=cell(1,length(UPstates));
for i=1:length(UPstates)
    currUPs=UPstates{i};
    currMUA=mua_ypoints{i};
    numUPs(i)=size(currUPs,1);
    thisTrialFRs=[];
    thisTrialdurs=[];
    for j=1:size(currUPs,1)
        if endsBy
            if currUPs(j,2)<=useUPsInTimeWindow(2)
            else
                continue
            end
        else
            if currUPs(j,1)>=useUPsInTimeWindow(1)
            else
                continue
            end
        end
        if (currUPs(j,2)-correctUPend)-currUPs(j,1)<UPdurationThresh
            continue
        end
        avFR=mean(currMUA(mua_xpoints>=currUPs(j,1) & mua_xpoints<=currUPs(j,2)));
        FRs=[FRs; avFR];
        thisTrialFRs=[thisTrialFRs; avFR];
        currDur=(currUPs(j,2)-correctUPend)-currUPs(j,1);
        durations=[durations; currDur];
        thisTrialdurs=[thisTrialdurs; currDur];
    end
    trialAvFR(i)=mean(thisTrialFRs);
    trialAvduration(i)=mean(thisTrialdurs);
    trialdurationHist{i}=thisTrialdurs;
end

bins=15;
[FR_heights,FR_centers]=histnorm(FRs,bins);
figure();
hist(FRs,bins);
bins=15;
[dur_heights,dur_centers]=histnorm(durations,bins);
figure(); 
hist(durations,bins);
bins=5;
[numUP_heights,numUP_centers]=histnorm(numUPs,bins);
figure(); 
hist(numUPs,bins);

FRhist.heights=FR_heights;
FRhist.centers=FR_centers;
FRhist.average=mean(FRs);
durhist.heights=dur_heights;
durhist.centers=dur_centers;
durhist.average=mean(durations);
numUPhist.heights=numUP_heights;
numUPhist.centers=numUP_centers;
numUPhist.average=mean(numUPs);

