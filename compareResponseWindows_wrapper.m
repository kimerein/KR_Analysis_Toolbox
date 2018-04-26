function runningResponseDiffs=compareResponseWindows_wrapper(LFPbySweep,ledForSweeps,LFP_Fs,ledAv)
% Wrapper for comparing response windows between LED and no LED conditions
%
% RETURNS:
% runningResponseDiffs: a vector of differences between the response window values
% for LEDdata and noLEDdata, taken as an average of multiple sub-samplings of the
% trial sets with and without LED pulses
% each row of responseDiffs represents a
% different response window, the mean value of noLEDdata for that response
% window is subtracted from the mean value of LEDdata for that response
% window
% 
% PARAMETERS:
% LFPbySweep: matrix with samples as columns and trials as rows
% responseWindows: an nx2 matrix specifying different response windows (in terms
% of time relative to the start of a trial), each row is a different
% response window; responseWindows(i,1) is the start of response window i,
% responseWindows(i,2) is the end of response window i
% ledForSweeps: a vector of 0's for no LED trials and 1's for LED trials,
% each element corresponds to a row of LFPbySweep
% LFP_Fs: sampling rate of LED and noLED data
% ledAv: a signal showing the LED for one (or averaged over all) trial(s),
% should be a vector with each element being a sample of the LED signal at
% a different time

baselineWindow=[0 0.2];
%baselineWindow=[0 0.1];
%baselineWindow=[0.5 2];
%baselineWindow=[0.5 1.89];
% %responseWindows=[2 2.05; 2.05 2.08; 2.08 2.1; 2.1 2.2; 2.2 2.4; 2.4 2.6; 2.6 2.8; 2.8 3; 3 3.2];
%responseWindows=[1.89 1.94; 1.94 1.97; 1.97 1.99; 1.99 2.09; 2.09 2.29; 2.29 2.49; 2.49 2.69; 2.69 2.89; 2.89 3.2];

%baselineWindow=[2.5 4];
%responseWindows=[2 2.05; 2.05 2.08; 2.08 2.1; 2.1 2.2; 2.2 2.4; 2.4 2.6; 2.6 2.8; 2.8 3; 3 3.2];
%responseWindows=[4 4.05; 4.05 4.08; 4.08 4.1; 4.1 4.2; 4.2 4.4; 4.4 4.6; 4.6 4.8; 4.8 5; 5 5.2];

responseWindows=[0.4 0.45;0.45 0.5;0.5 0.55;0.55 0.6;0.6 0.65;0.65 0.7;0.7 0.75;0.75 0.8;0.8 0.85;0.85 0.9;0.9 0.95;0.95 1;
                 1 1.05;1.05 1.1;1.1 1.15;1.15 1.2;1.2 1.25;1.25 1.3;1.3 1.35;1.35 1.4;1.4 1.45;1.45 1.5;1.5 1.55;1.55 1.6;1.6 1.65;1.65 1.7;1.7 1.75;
                 1.75 1.8;1.8 1.85;1.85 1.9;1.9 1.95;1.95 2;2 2.05;2.05 2.1;2.1 2.15;2.15 2.2;2.2 2.25;2.25 2.3;2.3 2.35;2.35 2.4;
                 2.4 2.45;2.45 2.5;2.5 2.55;2.55 2.6;2.6 2.65;2.65 2.7;2.7 2.75;2.75 2.8;2.8 2.85;2.85 2.9;2.9 2.95;2.95 3];

runningResponseDiffs=zeros(100,1+size(responseWindows,1));
for i=1:1
    responseDiffs=compareResponseWindows(LFPbySweep(ledForSweeps>0,:),LFPbySweep(~(ledForSweeps>0),:),responseWindows,baselineWindow,LFP_Fs,ledAv);
    runningResponseDiffs(i,:)=responseDiffs;
end
r=mean(runningResponseDiffs,1);
disp(r');
%disp(mean(runningResponseDiffs,1));
%disp(std(runningResponseDiffs,1));