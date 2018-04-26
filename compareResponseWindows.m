function [responseDiffs]=compareResponseWindows(LEDdata, noLEDdata, responseWindows, baselineWindow, LFP_Fs, ledAv)
% Note need to change save directory at bottom of this file, there are many
% parameters specified in this file, need to look at before using
% Calculates the differences in response windows between LED and no LED
% trials
% 
% RETURNS:
% responseDiffs: a vector of differences between the response window values
% for LEDdata and noLEDdata, each row of responseDiffs represents a
% different response window, the mean value of noLEDdata for that response
% window is subtracted from the mean value of LEDdata for that response
% window
% 
% PARAMETERS:
% LEDdata: matrix with samples as columns and trials as rows; only contains
% trials with LED on
% noLEDdata: matrix with samples as columns and trials as rows; only contains
% trials with LED off
% responseWindows: an nx2 matrix specifying different response windows (in terms
% of time relative to the start of a trial), each row is a different
% response window; responseWindows(i,1) is the start of response window i,
% responseWindows(i,2) is the end of response window i
% baselineWindow: a 2-element vector with the first element specifying the
% start of the time window to be used to calculate baseline, the second
% element specifying the end of the time window to be used to calculate
% baseline
% LFP_Fs: sampling rate of LED and noLED data
% ledAv: a signal showing the LED for one (or averaged over all) trial(s),
% should be a vector with each element being a sample of the LED signal at
% a different time

alignTo='baseline';
%alignTo='preDipPeriod';
saveSuffix='0';

% Match trial numbers
if size(LEDdata,1)>size(noLEDdata,1)
    takeInds=randperm(size(LEDdata,1));
    LEDdata=LEDdata(sort(takeInds(1:size(noLEDdata,1))),:);
elseif size(LEDdata,1)<size(noLEDdata,1)
    takeInds=randperm(size(noLEDdata,1));
    noLEDdata=noLEDdata(sort(takeInds(1:size(LEDdata,1))),:);
end

% Get baseline period values
baseInds=floor(baselineWindow(1)*LFP_Fs)+1:floor(baselineWindow(2)*LFP_Fs);

% Get averages
LEDav=mean(LEDdata,1);
noLEDav=mean(noLEDdata,1);

% Align averages to baseline
baseVal=mean(LEDav(baseInds));
LEDav=LEDav-baseVal;
baseVal=mean(noLEDav(baseInds));
noLEDav=noLEDav-baseVal;

% Align trials to baselines
% for i=1:size(LEDdata,1)
%     baseVal=mean(LEDdata(i,baseInds));
%     LEDdata(i,:)=LEDdata(i,:)-baseVal;
% end
% for i=1:size(noLEDdata,1)
%     baseVal=mean(noLEDdata(i,baseInds));
%     noLEDdata(i,:)=noLEDdata(i,:)-baseVal;
% end

LEDresponses=zeros(size(responseWindows,1),1);
noLEDresponses=zeros(size(responseWindows,1),1);
for i=1:size(responseWindows,1)
    inds=floor(responseWindows(i,1)*LFP_Fs)+1:floor(responseWindows(i,2)*LFP_Fs);
    LEDresponses(i)=mean(-LEDav(inds));
    noLEDresponses(i)=mean(-noLEDav(inds));
end

% Align responses to just prior to dip response
if ~strcmp(alignTo,'baseline')
    LEDbase=-LEDresponses(1);
    noLEDbase=-noLEDresponses(1);
    LEDresponses=LEDresponses-LEDresponses(1);
    noLEDresponses=noLEDresponses-noLEDresponses(1);
else
    LEDbase=0;
    noLEDbase=0.000001;
end
    
% Time vector
timeVec=0:responseWindows(end)/(size(LEDav,2)-1):responseWindows(end);

% Plot responses
responseXpoints=mean(responseWindows,2);
h=figure;
title('Responses for Different Response Windows');
subplot(4,1,1);
plot(timeVec,noLEDav,'Color','k');
hold on
plot(timeVec,LEDav,'Color','b');
min1=min(LEDav);
min2=min(noLEDav);
minall=min(min1,min2);
plot(timeVec,(ledAv/300)-min(ledAv/300)+minall,'Color','c');
axis([0 responseWindows(end) 0 1]);
axis tight
subplot(4,1,2);
plot([mean(baselineWindow);responseXpoints],[noLEDbase;noLEDresponses],'-','Color','k');
hold on
plot([mean(baselineWindow);responseXpoints],[LEDbase;LEDresponses],'--','Color','b');
axis([0 responseWindows(end) 0 1]);
axis 'auto y'
axis tight

h=figure;
title('LED Response as Fraction of No LED Response');
subplot(4,1,3);
fractionResponses=LEDresponses./noLEDresponses;
fractionResponses(1)=1;
plot([mean(baselineWindow);responseXpoints],[LEDbase/noLEDbase;fractionResponses],'-','Color','k');
axis([0 responseWindows(end) 0 1]);
axis 'auto y'
axis tight
saveas(h,['E:\Results\GammaLFPs\KR_2010-09-08\IntensityXDuration\LEDasFractionOfControlResponse' saveSuffix '.fig']);

figure;
subplot(4,1,4);
diffResponses=LEDresponses-noLEDresponses;
plot([mean(baselineWindow);responseXpoints],[LEDbase-noLEDbase;diffResponses],'-','Color','r');
%axis([0 responseWindows(end) -0.1 0.1]);
%axis 'auto y'
axis tight
disp('Average Response Difference');
disp(mean(diffResponses(2:end)));
saveas(h,['E:\Results\GammaLFPs\KR_2010-09-08\IntensityXDuration\diffResponseWindows' saveSuffix '.fig']);
saveLEDResponses=[LEDbase; LEDresponses];
savenoLEDResponses=[noLEDbase; noLEDresponses];
save(['E:\Results\GammaLFPs\KR_2010-09-08\IntensityXDuration\LEDResponses' saveSuffix '.mat'],'saveLEDResponses');
save(['E:\Results\GammaLFPs\KR_2010-09-08\IntensityXDuration\noLEDResponses' saveSuffix '.mat'],'savenoLEDResponses');

responseDiffs=[LEDbase-noLEDbase;diffResponses];