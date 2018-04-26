function [spikes,LFP,LFPPhoto_Fs,bandPassedLFP,photodiode,params]=alignToPhotodiode(spikes,LFP,LFPPhoto_Fs,bandPassedLFP,photodiode,params)
% Use of this function necessitates clean photodiode data
% LFP and photodiode data must have the same sampling rate
%
% Photo-aligns spikes, LFP and photodiode data
% This function aligns each trial of LFP data to the photodiode signal from
% that trial. Also, spike times are fixed to match the photo-aligned LFP.
% Columns of LFP are samples at different times relative to the
% onset of the trial, and rows of LFP are different trials. LFPPhoto_Fs is
% the sampling rate of both the photodiode data and LFP.
% bandPassedLFP is a matrix with the same set-up as LFP, however, the data
% represents a band-bassed version of LFP data. Also photo-align this
% band-passed data and returns in bandPassedLFP. photodiode contains the
% photodiode data for each trial. Rows are different trials; columns are
% difference times. params is a structure passed in containing
% information about the set-up of the trials. params must have a field
% params.ONstart, which contains the time, relative to the onset of a
% trial, when the stimulus turns on. After photo-aligning the data, this
% function will also fix params.ONstart to match the stimulus onset for the
% photo-aligned data.



% RIG-SPECIFIC!!!!!
% This is the time it takes the photodiode to charge from a min to a max
% = half a period of photodiode oscillation
chargingTime=13.75;   % in ms
chargingTime=chargingTime*10^-3;
chargingTimeInds=ceil(chargingTime*LFPPhoto_Fs);   % convert to indices  

% GET MOMENT OF STIMULUS ON SCREEN FROM PHOTODIODE DATA
% Absolute photodiode data is unimportant
% Make all photodiode data positive
% Smooth photodiode data
for i=1:size(photodiode,1)
    photodiode(i,:)=photodiode(i,:)-min(photodiode(i,:));
    for j=1:5
        photodiode(i,:)=smooth(photodiode(i,:)',3)';
    end
%   photodiode(i,:)=smooth(photodiode(i,:)',3)';

end

% During stimulus presentation, large photodiode current oscillations will
% be superimposed on high value
% Get local maxima and minima, take differences between adjacent max and
% min -> photodiode oscillations during stimulus will have greater
% amplitudes than some amplitude threshold
indsOfON=nan(size(photodiode,1),1);

for i=1:size(photodiode,1)
    if ~any(photodiode(i,:)~=0)
        figure;
        plot(photodiode(i,:));
        disp(['Photodiode trial is blank. Trial #: ' num2str(i)]);
        break
    end
    % Pre-process photodiode sweep
    % In order to use findpeaks function, need to eliminate steps or
    % multi-index peaks (identical neighboring values) in data
    [ppphotodiode,photoInds]=unique(photodiode(i,:),'last');
    [photoInds,sInds]=sort(photoInds);
    ppphotodiode=ppphotodiode(sInds);
    % Find threshold for separating high and low photodiode values
    % Note that TTL pulse always arrives before photodiode senses screen change
    % Thus, use initial photodiode value as a guess at its low value
    % Take first local max as photodiode high value
    % Get the low-high threshold as the average of these two values
    %lowValue=ppphotodiode(1);
    lowValue=min(ppphotodiode);
    %highValue=findpeaks(ppphotodiode,'npeaks',1);
    highValue=max(ppphotodiode);
    thresh=mean([lowValue highValue]);
    %thresh=photodiode(i,floor((3*10^-3)*LFPPhoto_Fs)+1);
%    thresh
    % Find last two peaks while photodiode active (stimulus on screen)
    [lastPks,lastPksInds]=findpeaks(fliplr(ppphotodiode),'npeaks',2,'minpeakheight',thresh);
    %[lastPks,lastPksInds]=findpeaks(fliplr(ppphotodiode),'npeaks',5,'minpeakheight',thresh);
    %lastPks=lastPks(1:2);
    %lastPksInds=lastPksInds(1:2);
%      figure;
%      plot(ppphotodiode);
%      hold on;
%      line([length(ppphotodiode)-fliplr(lastPksInds);
%      length(ppphotodiode)-fliplr(lastPksInds)],[1 1; 7 7]);
    % Compare min between these peaks to first neighboring peak
    % Difference between this min and this first peak used to calculate
    % threshold for photodiode oscillation amplitude when stimulus on
    % screen
    theseInds=(length(ppphotodiode)-(lastPksInds(2)-1)):(length(ppphotodiode)-(lastPksInds(1)-1));
    [interMins,interMinsInds]=findpeaks(max(ppphotodiode(theseInds))-ppphotodiode(theseInds),'npeaks',3,'sortstr','descend');
    interMin=ppphotodiode(theseInds(interMinsInds(1)));
    ampThresh=(1/2)*(lastPks(2)-interMin);
    peakDiffInds=theseInds(interMinsInds(1))-(length(ppphotodiode)-lastPksInds(2));
    %ampThresh
    % Get first peak with adjacent minima that are more than ampThresh 
    % current units away from the peak
    [allPks,allPksInds]=findpeaks(ppphotodiode,'minpeakheight',thresh);
    [allMins,allMinsInds]=findpeaks(max(ppphotodiode)-ppphotodiode);
    % Add first photodiode value as a minimum
    allMins=[ppphotodiode(1) ppphotodiode(allMinsInds)];
    allMinsInds=[1 allMinsInds];
    isMax=[ones(1,length(allPks)) zeros(1,length(allMins))];    % Combine local max and local min
    extInds=[allPksInds allMinsInds];
    extrema=[allPks allMins];
    [extInds,orderInds]=sort(extInds);
    isMax=isMax(orderInds);
    extrema=extrema(orderInds);
    % Find adjacent local min and max that are separated by more than
    % ampThresh
    firstPkInd=[];
    for j=find(isMax==1)
        if j==1
        elseif j+1>length(isMax) 
        else
            if isMax(j-1)==0 && isMax(j+1)==0
                %if extrema(j-1)+ampThresh<extrema(j) && extrema(j+1)+ampThresh<extrema(j)
                if extrema(j+1)+ampThresh<extrema(j)
                    if extrema(j-1)+ampThresh>=extrema(j)
                        firstPkInd=extInds(j+1)-peakDiffInds; % index into ppphotodiode
                        indsOfON(i)=photoInds(firstPkInd);
                    else
                        firstPkInd=extInds(j); % peak index into ppphotodiode
                        if photoInds(firstPkInd)-chargingTimeInds<1
                            indsOfON(i)=1;
                        else
                            % Convert back to photodiode indices from ppphotodiode indices
                            % and subtract photodiode charging time
                            indsOfON(i)=photoInds(firstPkInd)-chargingTimeInds;
                        end
                    end
                    break
%                     if extrema(j)<=peakNearVal
%                         firstPkInd=extInds(j); % peak index into ppphotodiode
%                     else
%                         firstPkInd=extInds(j);
%                         for k=extInds(j)+1:extInds(j+1) % find first time signal dips below peakNearVal  
%                             if ppphotodiode(k)<=peakNearVal
%                                 firstPkInd=k; % index into pphotodiode
%                                 break
%                             end
%                         end
%                     end
%                     break
                end
            end
        end
    end
%     figure;
%     plot(ppphotodiode);
%     hold on;
%     line([firstPkInd firstPkInd],[6 7],'Color','r');
    if isempty(firstPkInd)
        disp(['Warning: Could not find photodiode trigger for trial ' num2str(i)]);
        indsOfON(i)=NaN;
    end
        % Convert back to photodiode indices from ppphotodiode indices
        % and subtract photodiode charging time
%         if photoInds(firstPkInd)-chargingTimeInds<1
%             indsOfON(i)=1;
%         else
%             indsOfON(i)=photoInds(firstPkInd)-chargingTimeInds;
%         end
    %end
%     figure;
%     plot(photodiode(i,:));
%     hold on;
%     line([indsOfON(i)+chargingTimeInds indsOfON(i)+chargingTimeInds],[6 7],'Color','r');
%     'hi'
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
    offset=indsOfON(i)-startPointInd;
    if offset>0
        temp_LFP(:)=[LFP(i,offset+1:end) LFP(i,end)*ones(1,offset)];
        temp_bandPassedLFP(:)=[bandPassedLFP(i,offset+1:end) bandPassedLFP(i,end)*ones(1,offset)];
        temp_photodiode(:)=[photodiode(i,offset+1:end) photodiode(i,end)*ones(1,offset)];
    elseif offset<0
        temp_LFP(:)=[LFP(i,1)*ones(1,-offset) LFP(i,1:end+offset)];
        temp_bandPassedLFP(:)=[bandPassedLFP(i,1)*ones(1,-offset) bandPassedLFP(i,1:end+offset)];
        temp_photodiode(:)=[photodiode(i,1)*ones(1,-offset) photodiode(i,1:end+offset)];
    else
        temp_LFP(:)=LFP(i,:);
        temp_bandPassedLFP(:)=bandPassedLFP(i,:);
        temp_photodiode(:)=photodiode(i,:);
    end
    LFP(i,:)=temp_LFP;
    bandPassedLFP(i,:)=temp_bandPassedLFP;
    photodiode(i,:)=temp_photodiode;
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
figure;
plot((1:size(photodiode,2))'/LFPPhoto_Fs,photodiode);
m=max(photodiode(1,:));
hold on;
line([params.ONstart params.ONstart],[0 m],'Color','r');
title('Check alignment to photodiode');

