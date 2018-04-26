function [probability,duration,amplitude]=probabilityOfUpState(UPstates_LFP,trialLEDs,muax,muay)

LEDwindow=[1 2];
LEDwindow(2)=LEDwindow(2)-0.25;
LEDvals=[2 7 8 9];
noLEDvals=[1 3];

countWithLED=0;
countWithoutLED=0;
durationWithLED=nan(1,10*length(UPstates_LFP));
dLED_i=1;
durationWithoutLED=nan(1,10*length(UPstates_LFP));
dnoLED_i=1;
amplitudeWithLED=nan(1,10*length(UPstates_LFP));
aLED_i=1;
amplitudeWithoutLED=nan(1,10*length(UPstates_LFP));
anoLED_i=1;

for i=1:length(UPstates_LFP)
    curr=UPstates_LFP{i};
    for j=1:size(curr,1)
        if curr(1)>=LEDwindow(1) & curr(1)<=LEDwindow(2)
            % Begins during LED window
            % Is LED on?
            if ismember(trialLEDs(i),LEDvals)
                % Yes
                countWithLED=countWithLED+1;
                durationWithLED(dLED_i)=curr(2)-curr(1);
                dLED_i=dLED_i+1;
                m=muay{i};
                amplitudeWithLED(aLED_i)=nanmean(m(muax>=curr(1) & muax<=curr(2)));
                aLED_i=aLED_i+1;
            elseif ismember(trialLEDs(i),noLEDvals)
                % No
                countWithoutLED=countWithoutLED+1;
                durationWithoutLED(dnoLED_i)=curr(2)-curr(1);
                dnoLED_i=dnoLED_i+1;
                m=muay{i};
                amplitudeWithoutLED(anoLED_i)=nanmean(m(muax>=curr(1) & muax<=curr(2)));
                anoLED_i=anoLED_i+1;
            end
        end
    end
end

probabilityWithLED=countWithLED/sum(ismember(trialLEDs,LEDvals));
probabilityWithoutLED=countWithoutLED/sum(ismember(trialLEDs,noLEDvals));
probability.withLED=probabilityWithLED;
probability.withoutLED=probabilityWithoutLED;
probability.countWithLED=countWithLED;
probability.trialsWithLED=sum(ismember(trialLEDs,LEDvals));
probability.countWithoutLED=countWithoutLED;
probability.trialsWithoutLED=sum(ismember(trialLEDs,noLEDvals));

duration.withLED=durationWithLED;
duration.withoutLED=durationWithoutLED;

amplitude.withLED=amplitudeWithLED;
amplitude.withoutLED=amplitudeWithoutLED;

% figure(); 
% [n,xout]=hist(durationWithLED,10);
% plot(xout,n,'Color','k');
% hold on; 
% [n,xout]=hist(durationWithoutLED,10);
% plot(xout,n,'Color','r');
% title('Up state duration');

% figure(); 
% [n,xout]=hist(amplitudeWithLED,10);
% plot(xout,n,'Color','k');
% hold on; 
% [n,xout]=hist(amplitudeWithoutLED,10);
% plot(xout,n,'Color','r');
% title('Up state amplitude');
                