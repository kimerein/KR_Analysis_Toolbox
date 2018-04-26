function [A_vector,tau_vector,normOfResid_vector]=calculate_tau_across_layers(expt,useTheseTrodes,trodeChs,layerSpecs,waitS)

% useTheseTrodes should be a row vector of 4 trode indices
% these are the 4 trodes that contain the 16 ordered channels of 
% the linear probe
useTheseTrodes=[5];

% trodeChs should be a 4x4 matrix giving the 4 channels in each 
% trode of useTheseTrodes in each row
% trodeChs=[23    15    22    14;
%           19    11    20    12;
%           17     9    16     8;
%           18    10    21    13];
trodeChs=[23    15    22    14];

% layerSpecs should be a cell array of ints (1-16) specifying which
% channels in the 16 given by useTheseTrodes to include in each 
% cortical layer
% e.g.,
% layerSpecs{1}=[1 2]; layerSpecs{2}=[3 4 5 6]; etc.
layerSpecs{1}=[1 2 3];
% layerSpecs{2}=[4 5 6 7];
% layerSpecs{3}=[8 9];
% layerSpecs{4}=[10 11 12 13];
% layerSpecs{5}=[14 15 16];

% Different LED blocks and time of LED onset for each
% LEDblocks={3; [3 7]; 8};
% noLEDblocks={1; [1]; 1};
% LEDonsets={3.3; 3.3; 3.2};
LEDblocks={5};
noLEDblocks={0};
LEDonsets={1.2};

RigDef=RigDefs();
chSetup=RigDef.ChannelOrder{2}; % This is the 1x16 channel setup
for i=1:length(layerSpecs)
    layerSpecs{i}=chSetup(layerSpecs{i});
end

for i=1:length(layerSpecs)
    % Calculate taus for each cortical layer individually
    % First get spikes for this cortical layer
    theseChs=layerSpecs{i};
    currSpikes=[];
    for j=1:length(useTheseTrodes)
        if any(ismember(trodeChs(j,:),theseChs))
            evChsDefault=[1 2 3 4];
            includeChs=ismember(trodeChs(j,:),theseChs);
            includeEvChs=evChsDefault(includeChs);
            spikes=load([RigDef.Dir.Spikes expt.sort.trode(useTheseTrodes(j)).spikesfile]);
            spikes=spikes.spikes;
            spikes.event_channel=spikes.info.detect.event_channel';
            spikes=filtspikes(spikes,0,'fileInd',10:19);
            spikes=getTheseCh_spikes(spikes,includeEvChs);
            if isempty(currSpikes)
                currSpikes=spikes;
            else
                currSpikes=concatSpikes_noAssigns(currSpikes,spikes);
            end
        end
    end 
    spikes=currSpikes;
    block_A=zeros(length(LEDblocks),1);
    block_tau=zeros(length(LEDblocks),1);
    block_normOfResid=zeros(length(LEDblocks),1);
    for j=1:length(LEDblocks)
        % Now calculate tau using these spikes
        [xpoints1,ypoints1,ypoints2]=scriptForComparingMUA_addedParams(spikes,[],LEDblocks{j},noLEDblocks{j});
        [A,tau,normOfResid]=fitExponential(xpoints1,ypoints2,LEDonsets{j},waitS);
        block_A(j)=A;
        block_tau(j)=tau;
        block_normOfResid(j)=normOfResid;
    end
    A_vector{i}=block_A;
    tau_vector{i}=block_tau;
    normOfResid_vector{i}=block_normOfResid;
end


function theseChSpikes=getTheseCh_spikes(spikes,evChs)
% Get different channels
% spikes.event_channel=spikes.info.detect.event_channel';
theseChSpikes=filtspikes(spikes,0,'event_channel',evChs); 


function layerSpikes=getChbyCh_spikes(spikes)
% Get different channels
a=unique(spikes.info.detect.event_channel);
spikes.event_channel=spikes.info.detect.event_channel;
for i=1:length(a)
    layerSpikes{i}=filtspikes(spikes,0,'event_channel',a(i));
end


function [xpoints1,ypoints1,ypoints2]=getChbyCh_PSTHs(layerSpikes,ledValue,noLEDvalue)
for i=1:length(layerSpikes)
    [xpoints1{i},ypoints1{i},ypoints2{i}]=scriptForComparingMUA_addedParams(layerSpikes{i},[],ledValue,noLEDvalue);
end


    





% x_1=x1; y_1=y1; cx_1=x; cy_1=y;

% all_control_y=mean([cy_1;cy_2],1); all_expt_y=mean([y_1;y_2],1); all_control_x=cx_1; all_expt_x=x_1;

% [down_control_y,down_expt_y,down_x]=useLargerBinsize(all_control_y,all_expt_y,all_control_x,3);
% figure(); plot(all_control_x,all_control_y,'Color','k'); hold on; plot(all_expt_x,all_expt_y,'Color','r');
% [A,tau,normOfResid]=fitExponential(all_expt_x,all_expt_y,3.3);
% disp(tau);