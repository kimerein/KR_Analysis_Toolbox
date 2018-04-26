function [returnx_onset,returny_onset,returnx_led,returny_led,returnx_offset,returny_offset,usingTrials]=alignUPs_general(UPstates,trialLEDs,useLED,trialStims,useStimcond,LEDwindow,mua_x,mua_y)

% beforeLED=0.2; % in s
beforeLED=0.34; % in s
% beforeLED=0.2; % in s
bufferInds=5000;
takeBefore=0.2; % in s
takeAfter=1;
% takeBefore=0.5; % in s
% takeAfter=1;
% takeBefore=0.5; % in s
% takeAfter=0.5;
UPmax=10000;
showIndividualTrials=0;

if LEDwindow(1)-beforeLED-takeBefore<0
    disp('Switching takeBeforeOrAfter');
    takeBeforeOrAfter=LEDwindow(1)-beforeLED;
end

UPstates=UPstates(1:size(UPstates,1),1);

nIndsDuringLED=floor((LEDwindow(2)-LEDwindow(1))/(mua_x(2)-mua_x(1)));
alignStart=zeros(10000,bufferInds);
countStart=0;
alignLED=zeros(10000,bufferInds);
countLED=0;
LEDpdf=zeros(10000,bufferInds);
LEDpdf2=zeros(10000,bufferInds);
alignEnd=zeros(10000,bufferInds);
countEnd=0;
usingTrials=[];
for i=1:length(UPstates)
    currUPs=UPstates{i};
    trial_mua_y=mua_y{i};
    if ~ismember(trialLEDs(i),useLED) || ~ismember(trialStims(i),useStimcond)
        continue
    end
    for j=1:size(currUPs,1)
        % If this UP state begins within beforeLED of LED, then include
        % Another criterion for inclusion: UP state must continue at least
        % until LED starts
%         if currUPs(j,1)>=LEDwindow(1)-beforeLED && currUPs(j,2)>=LEDwindow(1)
%         if currUPs(j,1)>=LEDwindow(1)-beforeLED && currUPs(j,1)<=LEDwindow(1)
        if currUPs(j,1)>=LEDwindow(1)-beforeLED && currUPs(j,1)<=LEDwindow(1)
%             if max(trial_mua_y(mua_x>=LEDwindow(1)-beforeLED & mua_x<=LEDwindow(1)))>150 ...
%             || max(trial_mua_y(mua_x>=LEDwindow(1)-beforeLED & mua_x<=LEDwindow(1)))<50
%                 continue
%             end
            if currUPs(j,2)+takeAfter>max(mua_x)
                section_mua_y=trial_mua_y(mua_x>=currUPs(j,1)-takeBefore & mua_x<=max(mua_x));
                section_mua_x=mua_x(mua_x>=currUPs(j,1)-takeBefore & mua_x<=max(mua_x));
            else
                section_mua_y=trial_mua_y(mua_x>=currUPs(j,1)-takeBefore & mua_x<=currUPs(j,2)+takeAfter);
                section_mua_x=mua_x(mua_x>=currUPs(j,1)-takeBefore & mua_x<=currUPs(j,2)+takeAfter);
            end
            % Extra criterion -- mua_y can't get above UPmax
            if max(section_mua_y)>UPmax
                continue 
            end
%             if max(section_mua_y)>150
%                 continue 
%             end
            usingTrials=[usingTrials; i];
            % Align to start
            countStart=countStart+1;
            alignStart(countStart,1:length(section_mua_y))=section_mua_y;
            % Align to LED
            countLED=countLED+1;
            indAtLED=find(section_mua_x>=LEDwindow(1),1,'first');
            if isempty(indAtLED)
                indAtLED=find(mua_x>=LEDwindow(1),1,'first')-find(mua_x>=(currUPs(j,1)-takeBefore),1,'first');
            end
%             LEDpdf(countStart,indAtLED:indAtLED+nIndsDuringLED)=ones(size(indAtLED:indAtLED+nIndsDuringLED));
            LEDpdf(countStart,2500+indAtLED:2500+indAtLED+nIndsDuringLED)=ones(size(2500+indAtLED:2500+indAtLED+nIndsDuringLED));
            alignLED(countLED,(bufferInds/2)-indAtLED:(bufferInds/2)-indAtLED-1+length(section_mua_y))=section_mua_y;
            LEDpdf2(countLED,(bufferInds/2):(bufferInds/2)+nIndsDuringLED)=ones(size((bufferInds/2):(bufferInds/2)+nIndsDuringLED));
            % Align to end
            countEnd=countEnd+1;
            alignEnd(countEnd,size(alignEnd,2)-length(section_mua_y)+1:end)=section_mua_y;
        end
    end
end

timeStep=mua_x(2)-mua_x(1);
figure();
nInds=find(mean(alignStart(1:countStart-1,:),1)>0,1,'last');
nowx=0:timeStep:timeStep*(nInds-1);
plot(nowx,mean(alignStart(1:countStart-1,1:nInds),1));
hold on; 
ledpd=mean(LEDpdf(1:countStart-1,2500+1:2500+nInds),1);
plot(nowx,ledpd.*400,'Color','r');
title('UPs aligned to their onsets');
returnx_onset=nowx;
returny_onset=mean(alignStart(1:countStart-1,1:nInds),1);

figure(); 
firstInds=find(mean(alignLED(1:countLED-1,:),1)>0,1,'first');
lastInds=find(mean(alignLED(1:countLED-1,:),1)>0,1,'last');
nowx=firstInds*timeStep:timeStep:timeStep*(lastInds);
temp=nowx(1);
nowx=nowx-temp;
plot(nowx,mean(alignLED(1:countLED-1,firstInds:lastInds),1));
% plot(nowx,alignLED(1:countLED-1,firstInds:lastInds));
hold on; 
% line([(bufferInds/2)*timeStep-temp (bufferInds/2)*timeStep-temp],[0 500],'Color','r');
plot(nowx,mean(LEDpdf2(1:countLED-1,firstInds:lastInds),1).*400,'Color','r');
title('UPs aligned to the LED onset');
returnx_led=nowx;
returny_led=mean(alignLED(1:countLED-1,firstInds:lastInds),1);

if showIndividualTrials==1
    figure();
    plot(nowx,alignLED(1:countLED-1,firstInds:lastInds));
    hold on; 
    plot(nowx,mean(LEDpdf2(1:countLED-1,firstInds:lastInds),1).*400,'Color','r');
    title('UPs aligned to the LED onset -- individual trials');
end

figure(); 
nInds=find(mean(alignEnd(1:countEnd-1,:),1)>0,1,'first');
nowx=nInds*timeStep:timeStep:timeStep*(bufferInds);
nowx=nowx-nowx(1);
plot(nowx,mean(alignEnd(1:countEnd-1,nInds:end),1));
title('UPs aligned to their offsets');
returnx_offset=nowx;
returny_offset=mean(alignEnd(1:countEnd-1,nInds:end),1);