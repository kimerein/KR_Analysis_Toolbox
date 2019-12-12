% Script for analyzing FF data from Arbora

%% Add Chronux to path

addpath(genpath('chronux_2_11'));

%% Get daq file names

[theseAreDaqs,LED,Stim]=getDaqInDirectory('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-05\FF\DAQ files');

%% Get LFP

LFPchannel=2;
[LFPbySweep,Fs,completeSweeps]=getJustLFPbySweep('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-05\FF\DAQ files\',theseAreDaqs,20000,1,[LFPchannel]);
thetaDiff=getHippoTheta(LFPbySweep,zeros(1,size(LFPbySweep{1},1)),zeros(1,size(LFPbySweep{1},1)),0,0,20000);
noThetaTrials=nanmean(thetaDiff,2)<-0.03;

%% Get LED

% LEDchannel=22;
LEDchannel=23;
[LEDbySweep]=getJustLFPbySweep('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-05\FF\DAQ files\',theseAreDaqs,20000,1,[LEDchannel]);
temp=LEDbySweep{1};
LEDon=any(temp>0.5,2);

%% Load spikes and fix LED trials

load('C:\Users\Kim\Documents\MATLAB\from arbora gtacr2\kim sorted\AR2019-05-05_T3_Ch_21_20_19_18_spikes_SD4_kimsorted.mat');
trialsinexpt=unique(spikes.sweeps.trials);
ftrials=trialsinexpt(LEDon==1);

% ONLY USE FOR FF EXPERIMENTS
spikes.sweeps.led(ismember(spikes.sweeps.trials,ftrials))=spikes.sweeps.led(ismember(spikes.sweeps.trials,ftrials))+0.05;
spikes.led(ismember(spikes.trials,ftrials))=spikes.led(ismember(spikes.trials,ftrials))+0.05;
spikes.stimcond=spikes.led;
spikes.sweeps.stimcond=spikes.sweeps.led;

% find good units
temp=spikes.labels(:,1); 
useAssigns=temp(spikes.labels(:,2)==2); % 2 means good unit
disp(useAssigns);

%% Get PSTH for these units

trialDuration=4; % in seconds FOR FF
% trialDuration=7; % in seconds FOR DG
psth=measureAllUnitsResponseProperties(spikes,useAssigns,[0 trialDuration]);

%% Do F1 analysis

a=unique(spikes.trials);
uses=   {{[1 1.05]}; {[2 2.05]}; {[4 4.05]}; {[6 6.05]}; {[8 8.05]}; {[10 10.05]}; {[12 12.05]}; {[14 14.05]}; {[16 16.05]}; {[18 18.05]}; {[20 20.05]}; {[30 30.05]}; {[40 40.05]}; {[50 50.05]}; {[60 60.05]}};           
uses_tri={{a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}};

outputDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\AR2018-08-27 FF\T7';
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

doF1analysis_freqs([],[],outputDir,[],psth,freqs,freqs+0.05,[],uses,uses_tri,noThetaTrials);

%% Re-arrange by F1

arrangeF1byFreq([outputDir '\']);

%% Plot F1 analysis

plotThisFreq=6;
plotThisUnit=1;

plotDir=[outputDir '\Hz' num2str(plotThisFreq) '\'];

a=load([outputDir '\noTheta_LED.mat']);
times=a.noTheta.allS.t;

a=load([plotDir 'noTheta_trialAv_noLED.mat']);
F1_noTheta_noLED=a.noTheta_trialAv.F1amp;
a=load([plotDir 'noTheta_trialAv_LED.mat']);
F1_noTheta_LED=a.noTheta_trialAv.F1amp;
a=load([plotDir 'theta_trialAv_noLED.mat']);
F1_theta_noLED=a.theta_trialAv.F1amp;
a=load([plotDir 'theta_trialAv_LED.mat']);
F1_theta_LED=a.theta_trialAv.F1amp;

figure();
plot(times,F1_noTheta_noLED(plotThisUnit,:),'Color','k');
hold on;
plot(times,F1_noTheta_LED(plotThisUnit,:),'Color','c');
title('No Theta');

figure();
plot(times,F1_theta_noLED(plotThisUnit,:),'Color','k');
hold on;
plot(times,F1_theta_LED(plotThisUnit,:),'Color','c');
title('Theta');

plotWStderr_simple(F1_noTheta_noLED,F1_noTheta_LED,trialDuration,'k','c',times,1);
plotWStderr_simple(F1_theta_noLED,F1_theta_LED,trialDuration,'k','c',times,1);

%% Put together

plotThisFreq=6;
dd{1}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Reversing gratings\Mawake310 good thresh\F1 analysis' '\Hz' num2str(plotThisFreq)];
dd{2}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Reversing gratings\Mawake318 good thresh\F1 analysis' '\Hz' num2str(plotThisFreq)];
% dd{1}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Reversing gratings\Mawake311 good thresh\F1 analysis' '\Hz' num2str(plotThisFreq)];
% dd{2}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Reversing gratings\Mawake312 good thresh\F1 analysis' '\Hz' num2str(plotThisFreq)];
% dd{3}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Reversing gratings\Mawake313 good thresh\F1 analysis' '\Hz' num2str(plotThisFreq)];
% dd{4}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Reversing gratings\Mawake314 good thresh\F1 analysis' '\Hz' num2str(plotThisFreq)];
% dd{5}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Reversing gratings\Mawake315 good thresh\F1 analysis' '\Hz' num2str(plotThisFreq)];
% dd{6}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Reversing gratings\Mawake316 good thresh\F1 analysis' '\Hz' num2str(plotThisFreq)];

allLines=0;
plotF1analysis_simple(dd,trialDuration,times,allLines);
