function [spikes,LFP,LFPPhoto_Fs,bandPassedLFP,photodiode,params]=alignToPhotodiode(spikes,LFP,LFPPhoto_Fs,bandPassedLFP,photodiode,params)
% Use of this function necessitates clean photodiode data
% LFP and photodiode data must have the same sampling rate

% RIG-SPECIFIC!!!!!
% This is the time it takes the photodiode to charge from a min to a max
% = half a period of photodiode oscillation
chargingTime=11;   % in ms
chargingTimeInds=floor(chargingTime*LFPPhoto_Fs);   % convert to indices                

% GET MOMENT OF STIMULUS ON SCREEN FROM PHOTODIODE DATA
% Absolute photodiode data is unimportant
% Make all photodiode data positive
% Smooth photodiode data
for i=1:size(photodiode,1)
    photodiode(i,:)=photodiode(i,:)-min(photodiode(i,:));
    %photodiode(i,:)=smooth(photodiode(i,:)',10)';
end

% During stimulus presentation, large photodiode current oscillations will
% be superimposed on high value
% Get local maxima and minima, take differences between adjacent max and
% min -> photodiode oscillations during stimulus will have greater
% amplitudes than some amplitude threshold
indsOfON=zeros(size(photodiode,1),1);

for i=1:size(photodiode,1)
    % Pre-process photodiode sweep
    % In order to use findpeaks function, need to eliminate steps or
    % multi-index peaks (identical neighboring values) in data
    [ppphotodiode,photoInds]=unique(photodiode(i,:),'last');
    % Find threshold for separating high and low photodiode values
    % Note that TTL pulse always arrives before photodiode senses screen change
    % Thus, use initial photodiode value as a guess at its low value
    % Take first local max as photodiode high value
    % Get the low-high threshold as the average of these two values
    lowValue=ppphotodiode(1);
    highValue=findpeaks(ppphotodiode,'npeaks',1);
    thresh=mean([lowValue highValue]);
    % Find last two peaks while photodiode active (stimulus on screen)
    [lastPks,lastPksInds]=findpeaks(fliplr(ppphotodiode),'npeaks',2,'minpeakheight',thresh);
    % Compare min between these peaks to first neighboring peak
    % Difference between this min and this first peak used to calculate
    % threshold for photodiode oscillation amplitude when stimulus on
    % screen
    theseInds=(length(ppphotodiode)-(lastPks(2)-1)):(length(ppphotodiode)-(lastPks(1)-1));
    [interMins,interMinsInds]=findpeaks(max(ppphotodiode(theseInds))-ppphotodiode(theseInds),'npeaks',3,'sortstr','descend');
    interMin=ppphotodiode(interMinsInds(1));
    ampThresh=(1/3)*(lastPks(2)-interMin);
    % Get first peak with adjacent minima that are more than ampThresh 
    % current units away from the peak
    [allPks,allPksInds]=findpeaks(ppphotodiode,'minpeakheight',thresh);
    [allMins,allMinsInds]=findpeaks(max(ppphotodiode)-ppphotodiode);
    isMax=[ones(1,length(allPks)) zeros(1,length(allMins))];    % Combine local max and local min
    extInds=[allPksInds allMinsInds];
    [extInds,orderInds]=sort(extInds);
    isMax=isMax(orderInds);
    % Find adjacent local min and max that are separated by more than
    % ampThresh
    firstPkInd=[];
    for j=find(isMax==1)
        if j==1
        elseif j+1>length(isMax)
        else
            if isMax(j-1)==0 && isMax(j+1)==0
                if allMins(extInds(j-1))+ampThresh<allPks(extInds(j)) && allMins(extInds(j+1))+ampThresh<allPks(extInds(j))
                    firstPkInd=extInds(j); % peak index into ppphotodiode
                    break
                end
            end
        end
    end
    if isempty(firstPkInd)
        disp(['Warning: Could not find photodiode trigger for trial ' num2str(i));
        indsOfON(i)=NaN;
    else
        % Convert back to photodiode indices from ppphotodiode indices
        % and subtract photodiode charging time
        indsOfON(i)=photoInds(firstPkInd)-chargingTimeInds;
    end
    
    
    
    
    isMax=[ones(1,length(pks_locs)) zeros(1,length(mins_locs))];
    extremaVals=[pks photoMins];
    [combinedLocs,inds]=sort([pks_locs mins_locs]);
    isMax=isMax(inds);
    extremaVals=extremaVals(inds);
    % Find adjacent local min and max that are separated by more than
    % ampThresh; also, max must be greater than thresh
    bestMinInds=[];
    for j=find(isMax==0)
        if j+1>length(isMax)
            break
        end
        if isMax(j+1)==1
            if extremaVals(j+1)-extremaVals(j)>=ampThresh
                if extremaVals(j)~=lowValue
                    % Max-min difference big enough to be current
                    % oscillation during stim. & not initial charging to
                    % highValue
                    bestMinInds=[bestMinInds; combinedLocs(j)];
                end
            end
        end
    end
    %size(bestMinInds)
    dbstop error
    if ~isempty(bestMinInds)
        indsOfON(i)=min(bestMinInds);
    else
        indsOfON(i)=NaN;
    end
        
    for j=1:length(allPks)
        
    % Now, just subtract time for photodiode charging (should be constant across trials,
    % but is rig-specific)
    % PHOTODIODE CHARGING TIME IS RIG-SPECIFIC!!!!!!!!!!!!
    
    
    % If max value of photodiode for this sweep is NOT a local maximum,
    % suggests photodiode signal transiently saturated during TTL ON
    % In this case, find first "peak" of photodiode oscillation as the
    % last index with this max value
    sweepMax=max(ppphotodiode(:));
    [allPks,allPksInds]=findpeaks(photodiode(i,:),'minpeakheight',thresh,'sortstr','descend');
    if sweepMax~=
    
    
    % Get local maxima greater than highValue
    %plot(photodiode(i,:));
    %[pks,pks_locs]=findpeaks(photodiode(i,:),'threshold',highValue,'sortstr','descend');
%     [pks,pks_locs]=findpeaks(photodiode(i,:),'sortstr','descend');
% Get amplitude oscillation threshold from local minima closest to max
%    value of photodiode data
    %ampThresh=min(findpeaks(max(photodiode(i,:))-photodiode(i,:)));
    %disp('ampThresh')
    %ampThresh
    [pks,pks_locs]=findpeaks(photodiode(i,:),'sortstr','descend');  
    pks=pks(pks>thresh);
    pks_locs=pks_locs(pks>thresh);
    %disp('found peaks')
    %length(pks)
%     % Get local minima corresponding to local maxima
    [mins,mins_locs]=findpeaks(max(photodiode(i,:))-photodiode(i,:));
    photoMins=photodiode(i,mins_locs);
    % Combine local max and local min
    isMax=[ones(1,length(pks_locs)) zeros(1,length(mins_locs))];
    extremaVals=[pks photoMins];
    [combinedLocs,inds]=sort([pks_locs mins_locs]);
    isMax=isMax(inds);
    extremaVals=extremaVals(inds);
    % Find adjacent local min and max that are separated by more than
    % ampThresh; also, max must be greater than thresh
    bestMinInds=[];
    for j=find(isMax==0)
        if j+1>length(isMax)
            break
        end
        if isMax(j+1)==1
            if extremaVals(j+1)-extremaVals(j)>=ampThresh
                if extremaVals(j)~=lowValue
                    % Max-min difference big enough to be current
                    % oscillation during stim. & not initial charging to
                    % highValue
                    bestMinInds=[bestMinInds; combinedLocs(j)];
                end
            end
        end
    end
    %size(bestMinInds)
    dbstop error
    if ~isempty(bestMinInds)
        indsOfON(i)=min(bestMinInds);
    else
        indsOfON(i)=NaN;
    end
%    [pks,pks_locs]=findpeaks(photodiode(i,:),'threshold',ampThresh,'sortstr','descend');
    
end

% Align all indsOfON for LFP, bandPassedLFP, and photodiode
% Pad with zeros so don't have to worry about matching to totalTrialLength
% Very edges will not be correct in averages
startPointInd=floor(mean(indsOfON));
if startPointInd<1
    startPointInd=1;
end
for i=1:length(indsOfON)
    if isnan(indsOfON(i))
        disp(['WARNING: Could not find photodiode trigger for trial ' num2str(i)]);
        indsOfON(i)=startPointInd;
    end
end
for i=1:size(LFP,1)
    temp_LFP=zeros(1,size(LFP,2));
    temp_bandPassedLFP=zeros(1,size(bandPassedLFP,2));
    temp_photodiode=zeros(1,size(photodiode,2));
    offset=indsOfON(i)-startPointInd;
    if offset>0
        temp_LFP(:)=[LFP(i,offset+1:end) zeros(1,offset)];
        temp_bandPassedLFP(:)=[bandPassedLFP(i,offset+1:end) zeros(1,offset)];
        temp_photodiode(:)=[photodiode(i,offset+1:end) zeros(1,offset)];
    elseif offset<0
        temp_LFP(:)=[zeros(1,-offset) LFP(i,1:end-1)];
        temp_bandPassedLFP(:)=[zeros(1,-offset) bandPassedLFP(i,1:end-1)];
        temp_photodiode(:)=[zeros(1,-offset) photodiode(i,1:end-1)];
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
params.ONstart=startPointInd/LFPPhoto_Fs; 

% Plot aligned photodiode traces to check alignment
figure;
for i=1:size(photodiode,1)
    plot((1:size(photodiode,2))/LFPPhoto_Fs,photodiode(i,:));
    hold on;
    line([params.ONstart params.ONstart],[0 max(photodiode(1,:))],'Color','r');
end


