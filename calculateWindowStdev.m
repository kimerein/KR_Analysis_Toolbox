function calculateWindowStdev(LFPbySweep,ledForSweeps,LFP_Fs)
% Calculate the standard deviation of the response during each defined
% response window (see responseWindows below) for trials with and without
% an LED pulse
%
% PARAMETERS:
% LFPbySweep: matrix with samples as columns and trials as rows
% ledForSweeps: vector with 0's for trials without LED pulse, 1's for trials
% with LED pulse; rows correspond to trial rows of data
% LFP_Fs: sampling rate of LFPbySweep

LFPbySweep=LFPbySweep(1:end-1,:);
ledForSweeps=ledForSweeps(:,1:end-1);
% responseWindows=[0 0.2;0.4 0.45;0.45 0.5;0.5 0.55;0.55 0.6;0.6 0.65;0.65 0.7;0.7 0.75;0.75 0.8;0.8 0.85;0.85 0.9;0.9 0.95;0.95 1;
%                  1 1.05;1.05 1.1;1.1 1.15;1.15 1.2;1.2 1.25;1.25 1.3;1.3 1.35;1.35 1.4;1.4 1.45;1.45 1.5;1.5 1.55;1.55 1.6;1.6 1.65;1.65 1.7;1.7 1.75;
%                  1.75 1.8;1.8 1.85;1.85 1.9;1.9 1.95;1.95 2;2 2.05;2.05 2.1;2.1 2.15;2.15 2.2;2.2 2.25;2.25 2.3;2.3 2.35;2.35 2.4];

responseWindows=[0 0.2;0.4 0.45;0.45 0.5;0.5 0.55;0.55 0.6;0.6 0.65;0.65 0.7;0.7 0.75;0.75 0.8;0.8 0.85;0.85 0.9;0.9 0.95;0.95 1;
                 1 1.05;1.05 1.1;1.1 1.15;1.15 1.2;1.2 1.25;1.25 1.3;1.3 1.35;1.35 1.4;1.4 1.45;1.45 1.5;1.5 1.55;1.55 1.6;1.6 1.65;1.65 1.7;1.7 1.75;
                 1.75 1.8;1.8 1.85;1.85 1.9;1.9 1.95;1.95 2;2 2.05;2.05 2.1;2.1 2.15;2.15 2.2;2.2 2.25;2.25 2.3;2.3 2.35;2.35 2.4;
                 2.4 2.45;2.45 2.5;2.5 2.55;2.55 2.6;2.6 2.65;2.65 2.7;2.7 2.75;2.75 2.8;2.8 2.85;2.85 2.9;2.9 2.95;2.95 3];


LEDwindowStdevs=zeros(size(responseWindows,1),1);
for i=1:size(responseWindows,1)
    inds=floor(responseWindows(i,1)*LFP_Fs)+1:floor(responseWindows(i,2)*LFP_Fs);
    LEDwindowStdevs(i)=std(mean(LFPbySweep(ledForSweeps>0,inds),2));
end
noLEDwindowStdevs=zeros(size(responseWindows,1),1);
for i=1:size(responseWindows,1)
    inds=floor(responseWindows(i,1)*LFP_Fs)+1:floor(responseWindows(i,2)*LFP_Fs);
    noLEDwindowStdevs(i)=std(mean(LFPbySweep(~(ledForSweeps>0),inds),2));
end
disp(mean([noLEDwindowStdevs LEDwindowStdevs],2));