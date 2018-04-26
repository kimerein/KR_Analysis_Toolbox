function calculateLEDartifact_wrapper(daqDataName,daqFsName,ledForSweepsName)
% Wrapper for calculating the LED artifact
% Artifact is plotted and saved in calculateLEDartifact.m

% PARAMETERS:
% daqDataName: name of the .mat LFP data to use to calculate the artifact
% .mat file should contain a matrix with samples as columns and trials as
% rows
% daqFsName: name of the .mat sampling rate
% ledForSweepsName: name of the .mat file containing a vector with 1 for
% each trial with the LED on and 0 for each trial with the LED off

a=load(ledForSweepsName);
ledForSweeps=a.ledForSweeps;

AI=76;
includeTheseTrials=[41 55 59 63 67];
%includeTheseTrials=find(ledForSweeps>0);
LED_onset=1;
LED_offset=1.2;

a=load(daqDataName);
LFPbySweep=a.LFPbySweep;
a=load(daqFsName);
LFP_Fs=a.LFP_Fs;
calculateLEDartifact(LFPbySweep,LFP_Fs,LED_onset,LED_offset,includeTheseTrials,1,1,['E:\Results\GammaLFPs\KR_2010-08-23\LED Artifacts\AI' num2str(AI) '_LEDartifact.mat']);
