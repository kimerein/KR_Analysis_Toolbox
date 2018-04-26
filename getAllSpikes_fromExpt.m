function [concatSpikes, acrossTrodesChannelsNames]=getAllSpikes_fromExpt(expt, trode_names, fileInd)

% Note that this function is not symmetrical!
% It will only save in the concatenated structure the fields of the first
% encountered spikes struct

% Get and combine all spikes from all detected trodes in this experiment
% Return the subset of these spikes that occur in daq file numbers
% fileInd(1) through fileInd(2)
%
% expt          the expt struct specifying the experiment in which all these
%               units were recorded
% fileInd       the daq file indices to include in concatSpikes
%               if fileInd is [], all daq files will be included

RigDef=RigDefs();

concatSpikes={};
acrossTrodesChannels={};

thisExptTrodes={};
for i=1:length(expt.sort.trode)
    thisExptTrodes{i}=expt.sort.trode(i).name;
end
acrossTrodesSpikes={};
acrossTrodesChannelNames={};
z=1;
for i=1:length(trode_names)
    thisTrodeInd=0;
    for j=1:length(thisExptTrodes)
        if strcmp(trode_names(i),thisExptTrodes(j))
            thisTrodeInd=j;
        end
    end
    if thisTrodeInd==0
        disp('Error in matching selected trodes to trodes of this expt.');
        return
    end
    chsString='';
    for j=1:length(expt.sort.trode(thisTrodeInd).channels)
        chsString=strcat(chsString,num2str(expt.sort.trode(thisTrodeInd).channels(j)),'_');
    end
    spikesFile=strcat(RigDef.Dir.Spikes,expt.name,'_',trode_names(i),'_Ch_',chsString,'S1_spikes.mat');
    c=load(char(spikesFile));
    spikes=c.spikes;
    spikes.event_channel=spikes.info.detect.event_channel';
    for j=unique(spikes.info.detect.event_channel)'
        acrossTrodesChannelsNames{z}=expt.sort.trode(thisTrodeInd).channels(j);
        acrossTrodesSpikes{z}=filtspikes(filtspikes(spikes,0,'event_channel',j),0,'fileInd',fileInd);
        acrossTrodesSpikes{z}.event_channel=acrossTrodesChannelsNames{z}*ones(size(acrossTrodesSpikes{z}.event_channel));
        z=z+1;
    end      
end

fields=fieldnames(acrossTrodesSpikes{1});
concatSpikes=acrossTrodesSpikes{1};
for i=2:length(acrossTrodesSpikes)
    nfields=length(fieldnames(acrossTrodesSpikes{2}));
    if length(fields)~=nfields % should all have the same fields
        disp('The spikes structs you are trying to concatenate should all have the same fields!');
        concatSpikes=[];
        return
    end
    for j=1:length(fields)
        % If this field does not exist in one of the two spikes structures,
        % add it
        if ~isfield(acrossTrodesSpikes{i},fields{j})    % Don't care whether or not trode was sorted
            if strcmp(fields{j},'assigns')
                acrossTrodesSpikes{i}.assigns=nan(1,length(acrossTrodesSpikes{i}.spiketimes));
            elseif strcmp(fields{j},'labels')
                acrossTrodesSpikes{i}.labels=nan(1,2);
            end
        end
        if isstruct(acrossTrodesSpikes{i}.(fields{j}))
            insideFields=fieldnames(acrossTrodesSpikes{i}.(fields{j}));
            for k=1:length(insideFields)
%                 if isstruct(acrossTrodesSpikes{i}.(fields{j}).(insideFields{k}))
%                     inInsideFields=fieldnames(acrossTrodesSpikes{i}.(fields{j}).(insideFields{k}));
%                     for l=1:length(inInsideFields)
%                         doElse=1;
%                         if ~strcmp(insideFields{k},'display')
%                             if isa(concatSpikes.(fields{j}).(insideFields{k}).(inInsideFields{l}),'numeric')
%                                 concatSpikes.(fields{j}).(insideFields{k})=[concatSpikes.(fields{j}).(insideFields{k}).(inInsideFields{l}) acrossTrodesSpikes{i}.(fields{j}).(insideFields{k}).(inInsideFields{l})];
%                                 doElse=0;
%                             else
%                                 doElse=1;
%                             end
%                         end
%                         if doElse
%                             if i==2
%                                 temp=concatSpikes.(fields{j}).(insideFields{k}).(inInsideFields{l});
%                                 concatSpikes.(fields{j}).(insideFields{k}).(inInsideFields{l})={};
%                                 concatSpikes.(fields{j}).(insideFields{k}).(inInsideFields{l}){1}=temp;
%                                 concatSpikes.(fields{j}).(insideFields{k}).(inInsideFields{l}){2}=acrossTrodesSpikes{2}.(fields{j}).(insideFields{k}).(inInsideFields{l});
%                             else
%                                 concatSpikes.(fields{j}).(insideFields{k}).(inInsideFields{l}){i}=acrossTrodesSpikes{i}.(fields{j}).(insideFields{k}).(inInsideFields{l});
%                             end
%                         end
%                     end
                if isa(concatSpikes.(fields{j}).(insideFields{k}),'numeric')
                    concatSpikes.(fields{j}).(insideFields{k})=[concatSpikes.(fields{j}).(insideFields{k})  acrossTrodesSpikes{i}.(fields{j}).(insideFields{k})];
                else
%                 concatSpikes.(fields{j}).(insideFields{k})={concatSpikes.(fields{j}).(insideFields{k}) acrossTrodesSpikes{i}.(fields{j}).(insideFields{k})};
                    if i==2
                        temp=concatSpikes.(fields{j}).(insideFields{k});
                        concatSpikes.(fields{j}).(insideFields{k})={};
                        concatSpikes.(fields{j}).(insideFields{k}){1}=temp;
                        concatSpikes.(fields{j}).(insideFields{k}){2}=acrossTrodesSpikes{2}.(fields{j}).(insideFields{k});
                    else
                        concatSpikes.(fields{j}).(insideFields{k}){i}=acrossTrodesSpikes{i}.(fields{j}).(insideFields{k});
                    end
                end
            end
        else
            if size(acrossTrodesSpikes{i}.(fields{j}),1)>1
                % Pad the smaller waveforms with NaN
                if size(concatSpikes.(fields{j}),2)>size(acrossTrodesSpikes{i}.(fields{j}),2)
                    temp=[acrossTrodesSpikes{i}.(fields{j}) nan(size(acrossTrodesSpikes{i}.(fields{j}),1),size(concatSpikes.(fields{j}),2)-size(acrossTrodesSpikes{i}.(fields{j}),2),size(acrossTrodesSpikes{i}.(fields{j}),3))];
                    concatSpikes.(fields{j})=[concatSpikes.(fields{j}); temp];
                elseif size(concatSpikes.(fields{j}),2)<size(acrossTrodesSpikes{i}.(fields{j}),2)
                    temp=[concatSpikes.(fields{j}) nan(size(concatSpikes.(fields{j}),1),size(acrossTrodesSpikes{i}.(fields{j}),2)-size(concatSpikes.(fields{j}),2),size(concatSpikes.(fields{j}),3))];
                    concatSpikes.(fields{j})=[temp; acrossTrodesSpikes{i}.(fields{j})];
                else
                    concatSpikes.(fields{j})=[concatSpikes.(fields{j}); acrossTrodesSpikes{i}.(fields{j})];
                end
            else
                if strcmp(fields{j},'labels')
                    concatSpikes.(fields{j})=[concatSpikes.(fields{j}); acrossTrodesSpikes{i}.(fields{j})];
                else
                    if strcmp(fields{j},'event_channel')
                       concatSpikes.(fields{j})=[concatSpikes.(fields{j}) acrossTrodesSpikes{i}.(fields{j})];
                        %concatSpikes.(fields{j})=[concatSpikes.(fields{j}) i*acrossTrodesSpikes{i}.(fields{j})];
                    else
                        concatSpikes.(fields{j})=[concatSpikes.(fields{j}) acrossTrodesSpikes{i}.(fields{j})];
                    end
                end
            end
        end
    end
end
    
    
    
    
    