% Script for analyzing FF data from Arbora

%% Add Chronux to path

addpath(genpath('chronux_2_11'));

%% Get daq file names

[theseAreDaqs,LED,Stim]=getDaqInDirectory('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\arbora 2019-07-07\DG\DAQ');

% clear theseAreDaqs;
% for i=1:115
% theseAreDaqs{i}=['KR_2016-02-18_Mawake401_' num2str(i) '.daq'];
% end

%% Get LFP

LFPchannel=2; % arbora recordings
% LFPchannel=32; % kim recordings (config 2)
[LFPbySweep,Fs,completeSweeps]=getJustLFPbySweep('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\arbora 2019-07-07\DG\DAQ\',theseAreDaqs,20000,1,[LFPchannel]);
thetaDiff=getHippoTheta(LFPbySweep,zeros(1,size(LFPbySweep{1},1)),zeros(1,size(LFPbySweep{1},1)),0,0,20000);
noThetaTrials=nanmean(thetaDiff,2)<-0.05;
% noThetaTrials=nanmean(thetaDiff,2)<-0.1;
% noThetaTrials=nanmean(thetaDiff,2)<0.01376;

%% Get LED

LEDchannel=22;
% LEDchannel=6; % for some Kim expts
[LEDbySweep]=getJustLFPbySweep('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\arbora 2019-07-07\DG\DAQ\',theseAreDaqs,20000,1,[LEDchannel]);
temp=LEDbySweep{1};
LEDon=any(temp>0.5,2);

%% Load spikes and fix LED trials

load('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\arbora 2019-07-07\DG\Kim sorted\AR2019-07-07_T7_Ch_5_4_3_2_spikes_SD4.mat');
trialsinexpt=unique(spikes.sweeps.trials);
ftrials=trialsinexpt(LEDon==1);

% ONLY USE FOR FF EXPERIMENTS
% spikes.sweeps.led(ismember(spikes.sweeps.trials,ftrials))=spikes.sweeps.led(ismember(spikes.sweeps.trials,ftrials))+0.05;
% spikes.led(ismember(spikes.trials,ftrials))=spikes.led(ismember(spikes.trials,ftrials))+0.05;
% spikes.stimcond=spikes.led;
% spikes.sweeps.stimcond=spikes.sweeps.led;

% find good units
temp=spikes.labels(:,1); 
useAssigns=temp(spikes.labels(:,2)==2); % 2 means good unit
disp(useAssigns);

%% Get PSTH for these units

% trialDuration=4; % in seconds FOR FF
trialDuration=7; % in seconds FOR DG
% set LED conds in measureAllUnitsResponseProperties.m
psth=measureAllUnitsResponseProperties(spikes,useAssigns,[0 trialDuration]);

%% Do F1 analysis

a=unique(spikes.trials);
stico=unique(spikes.stimcond);
disp(stico);
uses={{[1]}; {[2]}; {[3]}; {[4]}; {[5]}; {[6]}; {[7]}; {[8]}};
uses_tri={{a}; {a}; {a}; {a}; {a}; {a}; {a}; {a}}; 
% uses={{[1]}; {[2]}; {[3]}; {[4]}};
% uses_tri={{a}; {a}; {a}; {a}}; 

outputDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\arbora 2019-07-07\DG\Sorted units F1 analysis\T7';

doF1analysis([],[],outputDir,[],psth,0,5,[],uses,uses_tri,noThetaTrials);


%% Plot F1 analysis

plotThisUnit=1;
temp=psth.psths{plotThisUnit}; 
figure(); 
plot(downSampAv(psth.t,10),downSampAv(nanmean(temp(psth.unitLED{plotThisUnit}==0,:),1),10),'Color','k'); 
hold on; 
plot(downSampAv(psth.t,10),downSampAv(nanmean(temp(psth.unitLED{plotThisUnit}==5,:),1),10),'Color','b');

plotDir=[outputDir '\'];

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

plotWStderr_simple(F1_noTheta_noLED,F1_noTheta_LED,trialDuration,'k','c',times,1); % all lines
plotWStderr_simple(F1_theta_noLED,F1_theta_LED,trialDuration,'k','c',times,1); % all lines

plotWStderr_simple(F1_noTheta_noLED,F1_noTheta_LED,trialDuration,'k','c',times,0); % average
plotWStderr_simple(F1_theta_noLED,F1_theta_LED,trialDuration,'k','c',times,0); % average

%% Put together

% plotThisFreq=4;
dd{1}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-05-05\Sorted units F1 analysis\T3_Ch_21_20_19_18';
dd{2}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-05-05\Sorted units F1 analysis\T4_Ch_17_16_15_14';
dd{3}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-05-05\Sorted units F1 analysis\T5_Ch_13_12_11_10';
dd{4}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-05-05\Sorted units F1 analysis\T6_Ch_9_8_7_6';
dd{5}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-05-05\Sorted units F1 analysis\T7_Ch_5_4_3_2';

dd{6}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-05\Sorted units F1 analysis\T4';
dd{7}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-05\Sorted units F1 analysis\T5';
dd{8}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-05\Sorted units F1 analysis\T6';
dd{9}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-05\Sorted units F1 analysis\T7';

dd{10}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-07\DG\Sorted units F1 analysis\T3';
dd{11}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-07\DG\Sorted units F1 analysis\T4';
dd{12}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-07\DG\Sorted units F1 analysis\T5';
dd{13}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-07\DG\Sorted units F1 analysis\T6';
dd{14}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-07-07\DG\Sorted units F1 analysis\T7';

dd{15}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-08-14\DG analyzed\Sorted units F1 analysis\T3';
dd{16}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-08-14\DG analyzed\Sorted units F1 analysis\T4';
dd{17}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-08-14\DG analyzed\Sorted units F1 analysis\T5';
dd{18}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-08-14\DG analyzed\Sorted units F1 analysis\T6';
dd{19}='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Gtacr2 final figure\AR2019-08-14\DG analyzed\Sorted units F1 analysis\T7';

allLines=0;
plotF1analysis_simple(dd,trialDuration,times,allLines);
