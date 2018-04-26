function test_UPs_dependence_on_thalamus(noLED_UPstates, LED_UPstates)

% times=0:0.1:6;
times=0:0.1:5;
ledWindow=[1.3 1.8]; % Time when LED was on during trial, in seconds

noLED_UPstates=noLED_UPstates(:,1);
LED_UPstates=LED_UPstates(:,1);

figure();
noledup=zeros(length(noLED_UPstates),length(times));
ledup=zeros(length(LED_UPstates),length(times));
for i=1:length(noLED_UPstates)
    currUPs=noLED_UPstates{i};
    for j=1:size(currUPs,1)
        noledup(i,times>currUPs(j,1) & times<currUPs(j,2))=1;
    end
end
for i=1:length(LED_UPstates)
    currUPs=LED_UPstates{i};
    for j=1:size(currUPs,1)
        ledup(i,times>currUPs(j,1) & times<currUPs(j,2))=1;
    end
end
plot(times,mean(noledup,1),'Color','k');
hold on; 
plot(times,mean(ledup,1),'Color','r');

% ledWindow=[1.05 1.55]; % Time when LED was on during trial, in seconds
% ledWindow=[3.6 4.1]; 
% ledWindow=[3.35 4.3]; 
% ledWindow=[1.05 2]; 
precedeTime=0.25; % in s
correctEndOfUPs=0.25;
continueDelay=0.1;

% What is the probability that an UP state begins during the LED?
% What is the probability that an UP state continues into the LED?
% UP states beginning during LED window
% Get stats:
% Count
% Duration
%
% UP states continuing into LED window
% Get stats:
% Count
% Duration
% Duration after LED onset
%
% UP states preceding LED window (ending <precedeTime before LED onset)
% Get stats:
% Count
% Interval before LED window
beginsInUP_numUPs=zeros(size(noLED_UPstates,1),1);
beginsInUP_durationUPs=[];
continuesIntoUP_numUPs=zeros(size(noLED_UPstates,1),1);
continuesIntoUP_durationUPs=[];
continuesIntoUP_durationAfterLED=[];
precedesUP_numUPs=zeros(size(noLED_UPstates,1),1);
precedesUP_interval=[];
for i=1:size(noLED_UPstates,1)
    currUPs=noLED_UPstates{i};
    for j=1:size(currUPs,1)
        if (currUPs(j,2)-correctEndOfUPs)-currUPs(j,1)<0.1 % Don't count UP states <100 ms in duration
            continue
        end
        % If this UP state ends during LED window
        if currUPs(j,2)-correctEndOfUPs>ledWindow(1) && currUPs(j,2)-correctEndOfUPs<=ledWindow(2)
            % If this UP state began during LED window
            if currUPs(j,1)>=ledWindow(1) && currUPs(j,1)<ledWindow(2)
                beginsInUP_numUPs(i)=beginsInUP_numUPs(i)+1;
                beginsInUP_durationUPs=[beginsInUP_durationUPs; (currUPs(j,2)-correctEndOfUPs)-currUPs(j,1)];
            else
                % means this UP state began before LED window
                continuesIntoUP_numUPs(i)=continuesIntoUP_numUPs(i)+1;
                continuesIntoUP_durationUPs=[continuesIntoUP_durationUPs; (currUPs(j,2)-correctEndOfUPs)-currUPs(j,1)];
                continuesIntoUP_durationAfterLED=[continuesIntoUP_durationAfterLED; (currUPs(j,2)-correctEndOfUPs)-ledWindow(1)];
            end
        else
            % If UP state began before LED window
            if currUPs(j,1)<ledWindow(1)
                % If UP state ends after LED window
                if currUPs(j,2)-correctEndOfUPs>ledWindow(2)
                    continuesIntoUP_numUPs(i)=continuesIntoUP_numUPs(i)+1;
                    continuesIntoUP_durationUPs=[continuesIntoUP_durationUPs; (currUPs(j,2)-correctEndOfUPs)-currUPs(j,1)];
                    continuesIntoUP_durationAfterLED=[continuesIntoUP_durationAfterLED; (currUPs(j,2)-correctEndOfUPs)-ledWindow(1)];
                else
                    % means UP state ends before LED window
                    % If UP state ends within precedeTime of onset of LED window
                    if ledWindow(1)-(currUPs(j,2)-correctEndOfUPs)<=precedeTime
                        precedesUP_numUPs(i)=precedesUP_numUPs(i)+1;
                        precedesUP_interval=[precedesUP_interval; ledWindow(1)-(currUPs(j,2)-correctEndOfUPs)];
                    else
                        % means UP state ends before LED window
                    end
                end
            elseif currUPs(j,1)>=ledWindow(1) && currUPs(j,1)<ledWindow(2)
                % means UP state began during LED window but ends after it
                disp('currUPs begins during LED and ends after LED');
                    disp(currUPs(j,:));
                beginsInUP_numUPs(i)=beginsInUP_numUPs(i)+1;
                beginsInUP_durationUPs=[beginsInUP_durationUPs; (currUPs(j,2)-correctEndOfUPs)-currUPs(j,1)];
            end
        end
    end 
end
noLED.beginsInUP_numUPs=beginsInUP_numUPs;
noLED.beginsInUP_durationUPs=beginsInUP_durationUPs;
noLED.continuesIntoUP_numUPs=continuesIntoUP_numUPs;
noLED.continuesIntoUP_durationUPs=continuesIntoUP_durationUPs;
noLED.continuesIntoUP_durationAfterLED=continuesIntoUP_durationAfterLED;
noLED.precedesUP_numUPs=precedesUP_numUPs;
noLED.precedesUP_interval=precedesUP_interval;
disp('Fraction of trials without LED in which UP state begins during ledWindow');
disp(sum(noLED.beginsInUP_numUPs>0 & noLED.precedesUP_numUPs<0.5));
disp('over');
disp(length(noLED.beginsInUP_numUPs));
disp('equals');
disp(sum(noLED.beginsInUP_numUPs>0 & noLED.precedesUP_numUPs<0.5)/length(noLED.beginsInUP_numUPs));
disp('Given that UP state continues into ledWindow, probability that it continues for >100 ms after LED onset');
disp(sum(noLED.continuesIntoUP_durationAfterLED>0.1));
disp('over');
disp(length(noLED.continuesIntoUP_durationAfterLED));
disp('equals');
disp(sum(noLED.continuesIntoUP_durationAfterLED>0.1)/length(noLED.continuesIntoUP_durationAfterLED));
                  
beginsInUP_numUPs=zeros(size(LED_UPstates,1),1);
beginsInUP_durationUPs=[];
continuesIntoUP_numUPs=zeros(size(LED_UPstates,1),1);
continuesIntoUP_durationUPs=[];
continuesIntoUP_durationAfterLED=[];
precedesUP_numUPs=zeros(size(LED_UPstates,1),1);
precedesUP_interval=[];
for i=1:size(LED_UPstates,1)
    currUPs=LED_UPstates{i};
    for j=1:size(currUPs,1)
        % If this UP state ends during LED window
        if currUPs(j,2)>ledWindow(1) && currUPs(j,2)<=ledWindow(2)
            % If this UP state began during LED window
            if currUPs(j,1)>=ledWindow(1) && currUPs(j,1)<ledWindow(2)
                beginsInUP_numUPs(i)=beginsInUP_numUPs(i)+1;
                beginsInUP_durationUPs=[beginsInUP_durationUPs; currUPs(j,2)-currUPs(j,1)];
            else
                % means this UP state began before LED window
                continuesIntoUP_numUPs(i)=continuesIntoUP_numUPs(i)+1;
                continuesIntoUP_durationUPs=[continuesIntoUP_durationUPs; currUPs(j,2)-currUPs(j,1)];
                continuesIntoUP_durationAfterLED=[continuesIntoUP_durationAfterLED; currUPs(j,2)-ledWindow(1)];
            end
        else
            % If UP state began before LED window
            if currUPs(j,1)<ledWindow(1)
                % If UP state ends after LED window
                if currUPs(j,2)>ledWindow(2)
                    continuesIntoUP_numUPs(i)=continuesIntoUP_numUPs(i)+1;
                    continuesIntoUP_durationUPs=[continuesIntoUP_durationUPs; currUPs(j,2)-currUPs(j,1)];
                    continuesIntoUP_durationAfterLED=[continuesIntoUP_durationAfterLED; currUPs(j,2)-ledWindow(1)];
                else
                    % means UP state ends before LED window
                    % If UP state ends within precedeTime of onset of LED window
                    if ledWindow(1)-currUPs(j,2)<=precedeTime
                        precedesUP_numUPs(i)=precedesUP_numUPs(i)+1;
                        precedesUP_interval=[precedesUP_interval; ledWindow(1)-currUPs(j,2)];
                    else
                        % means UP state ends before LED window
                    end
                end
            elseif currUPs(j,1)>=ledWindow(1) && currUPs(j,1)<ledWindow(2)
                % means UP state began during LED window but ends after it
                beginsInUP_numUPs(i)=beginsInUP_numUPs(i)+1;
                beginsInUP_durationUPs=[beginsInUP_durationUPs; currUPs(j,2)-currUPs(j,1)];
            end
        end
    end 
end
LED.beginsInUP_numUPs=beginsInUP_numUPs;
LED.beginsInUP_durationUPs=beginsInUP_durationUPs;
LED.continuesIntoUP_numUPs=continuesIntoUP_numUPs;
LED.continuesIntoUP_durationUPs=continuesIntoUP_durationUPs;
LED.continuesIntoUP_durationAfterLED=continuesIntoUP_durationAfterLED;
LED.precedesUP_numUPs=precedesUP_numUPs;
LED.precedesUP_interval=precedesUP_interval;
disp('Fraction of trials WITH LED in which UP state begins during ledWindow');
disp(sum(LED.beginsInUP_numUPs>0 & LED.precedesUP_numUPs<0.5));
disp('over');
disp(length(LED.beginsInUP_numUPs));
disp('equals');
disp(sum(LED.beginsInUP_numUPs>0 & LED.precedesUP_numUPs<0.5)/length(LED.beginsInUP_numUPs));
disp('Given that UP state continues into ledWindow, probability that it continues for >100 ms after LED onset');
disp(sum(LED.continuesIntoUP_durationAfterLED>0.1));
disp('over');
disp(length(LED.continuesIntoUP_durationAfterLED));
disp('equals');
disp(sum(LED.continuesIntoUP_durationAfterLED>0.1)/length(LED.continuesIntoUP_durationAfterLED));

% Make figures comparing ledWindow with and without LED
% beginsInUP_numUPS
bins=10;
makeComparisonFigure(noLED.beginsInUP_numUPs,LED.beginsInUP_numUPs,bins);
title('beginsInUP_numUPs');
bins=10;
makeComparisonFigure(noLED.beginsInUP_durationUPs,LED.beginsInUP_durationUPs,bins);
title('beginsInUP_durationUPs');
bins=10;
makeComparisonFigure(noLED.continuesIntoUP_numUPs,LED.continuesIntoUP_numUPs,bins);
title('continuesIntoUP_numUPs');
bins=10;
makeComparisonFigure(noLED.continuesIntoUP_durationUPs,LED.continuesIntoUP_durationUPs,bins);
title('continuesIntoUP_durationUPs');
bins=10;
makeComparisonFigure(noLED.continuesIntoUP_durationAfterLED,LED.continuesIntoUP_durationAfterLED,bins);
title('continuesIntoUP_durationAfterLED');
% bins=10;
% makeComparisonFigure(noLED.precedesUP_numUPs,LED.precedesUP_numUPs,bins);
% title('precedesUP_numUPs');
% bins=10;
% makeComparisonFigure(noLED.precedesUP_interval,LED.precedesUP_interval,bins);
% title('precedesUP_interval');
end

function makeComparisonFigure(noLED,LED,bins)

[heights,centers]=histnorm(noLED,bins);
figure(); 
plot(centers,heights,'Color','k');
hold on;
[heights,centers]=histnorm(LED,bins);
plot(centers,heights,'Color','r');

end

            