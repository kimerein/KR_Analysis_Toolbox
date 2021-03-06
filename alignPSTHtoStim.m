function [newPSTH,ledConds,LEDmode]=alignPSTHtoStim(psth,ledLFP,ledConds,ledVal,times,LEDmode)

newPSTH=[];
% ledThresh=2.5;
ledThresh=0.2;
% ledThresh=1.5;
alignLEDoffTrialsToo=false;

if iscell(ledLFP)
    l=ledLFP{1};
else
    l=ledLFP;
end

if any(isnan(ledConds))
    l=l(~isnan(ledConds),:);
    ledConds=ledConds(~isnan(ledConds));
end

if length(ledConds)~=size(psth.psths{1},1)
    disp('lengths of ledConds and psth.psths do not match!');
    return
end

LEDon=l(ismember(single(ledConds),single(ledVal)),:)>ledThresh;
LEDstarts=zeros(sum(ismember(single(ledConds),single(ledVal))),1);
ledOnTrials=find(ismember(single(ledConds),single(ledVal)));
for i=1:size(LEDon,1)
    temp=find((LEDon(i,:)==1),1,'first');
    if isempty(temp)
        ledConds(ledOnTrials(i))=0;
        LEDstarts(i)=nan;
    else
        LEDstarts(i)=temp;
    end
end
LEDstarts=LEDstarts(~isnan(LEDstarts));
LEDstartsTime=zeros(length(LEDstarts),1);
for i=1:length(LEDstarts)
    LEDstartsTime(i)=times(LEDstarts(i));
end

if isempty(LEDmode)
    LEDmode=mode(LEDstartsTime);
end

% Align PSTH
newPSTH.t=psth.t;
newPSTH.psths=cell(length(psth.psths),1);
newPSTH.unitTrials=psth.unitTrials;
newPSTH.unitStimcond=psth.unitStimcond;
newPSTH.unitLED=psth.unitLED;
timeInd=psth.t(2)-psth.t(1);
checkLED=l;
n=10;
ledOnTrials=find(ismember(single(ledConds),single(ledVal)));
ledOffTrials=find(~ismember(single(ledConds),single(ledVal)));
timeIndRatio=(timeInd/n)./(times(2)-times(1));
for j=1:length(newPSTH.unitLED)
    currunitled=newPSTH.unitLED{j};
    currunitled(~ismember(1:length(currunitled),ledOnTrials))=0;
    newPSTH.unitLED{j}=currunitled;
end
for j=1:length(psth.psths)
    currp=psth.psths{j};
    for i=1:length(LEDstartsTime)
        currTime=LEDstartsTime(i);
        currTrial=ledOnTrials(i);
        if currTime~=LEDmode
            shiftby=LEDmode-currTime;
            shiftbybackup=shiftby;
            shiftby=abs(shiftby);
            inInds=shiftby/timeInd;
            remainder=inInds-floor(inInds);
            rem=round(remainder*10);
            rem=rem*(n/10);
            inInds=floor(inInds);
            interp_vec=interp(currp(currTrial,:),n);
            if shiftbybackup>0 && inInds>0
                shifted_interp_vec=[zeros(1,inInds*n+rem) interp_vec(1:end-(inInds*n+rem))];
                if j==1
                    checkLED(currTrial,:)=[zeros(1,floor((inInds*n+rem)*timeIndRatio)) checkLED(currTrial,1:end-floor((inInds*n+rem)*timeIndRatio))];
                end
            elseif inInds>0
                shifted_interp_vec=[interp_vec(inInds*n+rem+1:end) zeros(1,inInds*n+rem)];
                if j==1
                    checkLED(currTrial,:)=[checkLED(currTrial,floor((inInds*n+rem)*timeIndRatio)+1:end) zeros(1,floor((inInds*n+rem)*timeIndRatio))];
                end
            else
                shifted_interp_vec=interp_vec;
            end
            ds_shifted_interp_vec=shifted_interp_vec(1:n:end);
            currp(currTrial,:)=ds_shifted_interp_vec;
        end
    end
    if alignLEDoffTrialsToo==true
        LEDOffstartsTime=nan(1,length(ledOffTrials));
        for i=1:length(ledOffTrials)
            [~,mi]=nanmin(abs(ledOffTrials(i)-ledOnTrials)); % find closest led on trial
            LEDOffstartsTime(i)=LEDstartsTime(mi); % shift led off trial by same amount
        end
        for i=1:length(ledOffTrials)
            currTime=LEDOffstartsTime(i);
            currTrial=ledOffTrials(i);
            if currTime~=LEDmode
                shiftby=LEDmode-currTime;
                shiftbybackup=shiftby;
                shiftby=abs(shiftby);
                inInds=shiftby/timeInd;
                remainder=inInds-floor(inInds);
                rem=round(remainder*10);
                rem=rem*(n/10);
                inInds=floor(inInds);
                interp_vec=interp(currp(currTrial,:),n);
                if shiftbybackup>0 && inInds>0
                    shifted_interp_vec=[zeros(1,inInds*n+rem) interp_vec(1:end-(inInds*n+rem))];
                    if j==1
                        checkLED(currTrial,:)=[zeros(1,floor((inInds*n+rem)*timeIndRatio)) checkLED(currTrial,1:end-floor((inInds*n+rem)*timeIndRatio))];
                    end
                elseif inInds>0
                    shifted_interp_vec=[interp_vec(inInds*n+rem+1:end) zeros(1,inInds*n+rem)];
                    if j==1
                        checkLED(currTrial,:)=[checkLED(currTrial,floor((inInds*n+rem)*timeIndRatio)+1:end) zeros(1,floor((inInds*n+rem)*timeIndRatio))];
                    end
                else
                    shifted_interp_vec=interp_vec;
                end
                ds_shifted_interp_vec=shifted_interp_vec(1:n:end);
                currp(currTrial,:)=ds_shifted_interp_vec;
            end
        end
    end
    newPSTH.psths{j}=currp;
end

figure();
plot(downSampAv(times,10),downSampMatrix(l,10)');
title('Before');

figure();
plot(downSampAv(times,10),downSampMatrix(checkLED,10)');
title('After');

% Throw out non-aligned trials


end
    