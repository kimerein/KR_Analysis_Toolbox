function [spikes,LFP,LFPPhoto_Fs,bandPassedLFP,photodiode,ledBySweep,params]=photoAlignData(spikes,LFP,LFPPhoto_Fs,bandPassedLFP,photodiode,ledBySweep,params)
% Use of this function necessitates clean photodiode data
% LFP and photodiode data must have the same sampling rate

global saveToDir2

% RIG-SPECIFIC!!!!!
% Config. for Kim's rig:
% Threshold for photodiode ON detection
% When photodiode signal dips below this value, photodiode ON, align to
% this X point
%y_thresh=-0.52;
%y_thresh=-0.23;

% Find times when photodiode turns ON for each trial
indsOfON=nan(size(photodiode,1),1);
for i=1:size(photodiode,1)
    %if photodiode(i,1)>-0.3 % For gain=500
        %y_thresh=-0.23;
        %y_thresh=-0.60;
        y_thresh=-0.65;
        %y_high_thresh=-0.46;
        %y_high_thresh=-0.58;
        y_high_thresh=-0.6;
%     if photodiode(i,1)>-0.35 % For gain=200
%         y_thresh=-0.22;
%         y_high_thresh=-0.46;
    %else
    %    y_thresh=-0.52;
    %    y_high_thresh=-0.46;
    %end
    if ~any(photodiode(i,:)~=0)
        figure;
        plot(photodiode(i,:));
        disp(['Photodiode trial is blank. Trial #: ' num2str(i)]);
        break
    end
   % In case sweep starts with signal value low, find first data point with
   % y-value greater than y_thresh, then find the next data point with y-value 
   % less than y_thresh
   first_high_ind=find(photodiode(i,:)>y_high_thresh,1,'first');
   first_below_ind=find(photodiode(i,first_high_ind:end)<=y_thresh,1,'first')+first_high_ind-1;
   if isempty(first_below_ind)
       disp(['Warning: Trial # ' num2str(i) ' does not have photodiode trigger.']);
   else
       indsOfON(i)=first_below_ind;
   end
end

% figure;
% plot(photodiode');
% hold on;
% line([indsOfON'+chargingTimeInds; indsOfON'+chargingTimeInds],[ones(1,length(indsOfON))*6; ones(1,length(indsOfON))*7]);

% figure;
% plot((1:size(photodiode,2))/LFPPhoto_Fs,photodiode);
% hold on;
% line([indsOfON';indsOfON']/LFPPhoto_Fs,[ones(1,length(indsOfON))*6; ones(1,length(indsOfON))*7]);

% Align all indsOfON for LFP, bandPassedLFP, and photodiode
% Pad with nearest LFP value so don't have to worry about matching to totalTrialLength
% Very edges will not be correct in averages
startPointInd=floor(mean(indsOfON));
if startPointInd<1
    startPointInd=1;
end
indsOfON(isnan(indsOfON))=startPointInd*ones(sum(isnan(indsOfON)),1);
for i=1:size(LFP,1)
    temp_LFP=zeros(1,size(LFP,2));
    temp_bandPassedLFP=zeros(1,size(bandPassedLFP,2));
    temp_photodiode=zeros(1,size(photodiode,2));
    temp_ledBySweep=zeros(1,size(photodiode,2));
    offset=indsOfON(i)-startPointInd;
    if offset>0
        temp_LFP(:)=[LFP(i,offset+1:end) LFP(i,end)*ones(1,offset)];
        temp_bandPassedLFP(:)=[bandPassedLFP(i,offset+1:end) bandPassedLFP(i,end)*ones(1,offset)];
        temp_photodiode(:)=[photodiode(i,offset+1:end) photodiode(i,end)*ones(1,offset)];
        temp_ledBySweep(:)=[ledBySweep(i,offset+1:end) ledBySweep(i,end)*ones(1,offset)];
    elseif offset<0
        temp_LFP(:)=[LFP(i,1)*ones(1,-offset) LFP(i,1:end+offset)];
        temp_bandPassedLFP(:)=[bandPassedLFP(i,1)*ones(1,-offset) bandPassedLFP(i,1:end+offset)];
        temp_photodiode(:)=[photodiode(i,1)*ones(1,-offset) photodiode(i,1:end+offset)];
        temp_ledBySweep(:)=[ledBySweep(i,1)*ones(1,-offset) ledBySweep(i,1:end+offset)];
    else
        temp_LFP(:)=LFP(i,:);
        temp_bandPassedLFP(:)=bandPassedLFP(i,:);
        temp_photodiode(:)=photodiode(i,:);
        temp_ledBySweep(:)=ledBySweep(i,:);
    end
    LFP(i,:)=temp_LFP;
    bandPassedLFP(i,:)=temp_bandPassedLFP;
    photodiode(i,:)=temp_photodiode;
    ledBySweep(i,:)=temp_ledBySweep;
end
    
% Also adjust spike times to match aligned stimulus onsets
for i=unique(spikes.trials)
    timeOffset=(indsOfON(i)-startPointInd)/LFPPhoto_Fs;
    if offset>0
        % Subtract timeOffset from all spiketimes for trial i
        spikes.spiketimes(spikes.trials==i)=spikes.spiketimes(spikes.trials==i)-timeOffset;
        % Any spikes at time < 0 -> just put at time = 0
        spikes.spiketimes(spikes.spiketimes<0)=0;
    elseif offset<0
        % Add timeOffset to all spiketimes for trial i
        spikes.spiketimes(spikes.trials==i)=spikes.spiketimes(spikes.trials==i)+timeOffset;
        % Any spikes at time > params.totalTrialLength -> just put at time
        % params.totalTrialLength
        spikes.spiketimes(spikes.spiketimes>params.totalTrialLength)=params.totalTrialLength;
    end
end

% Adjust stimulus parameters
params.ONstart=(startPointInd-1)/LFPPhoto_Fs; 

% Plot aligned photodiode traces to check alignment
thisFig=figure;
plot((1:size(photodiode,2))'/LFPPhoto_Fs,photodiode);
m=max(photodiode(1,:));
hold on;
line([params.ONstart params.ONstart],[0 m],'Color','r');
title('Check alignment to photodiode');
saveas(thisFig,strcat(saveToDir2,'\checkPhotoalign.fig'));

