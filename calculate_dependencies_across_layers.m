function [layers]=calculate_dependencies_across_layers(expt,useTheseTrodes,trodeChs,layerSpecs,F1amps)

% useTheseTrodes should be a row vector of 4 trode indices
% these are the 4 trodes that contain the 16 ordered channels of 
% the linear probe
useTheseTrodes=[1 2 3 4];

% use these file inds
useFileInds=[55:83];

% use these stim conds
% useStimConds=1:2;

% trodeChs should be a 4x4 matrix giving the 4 channels in each 
% trode of useTheseTrodes in each row
trodeChs=[23    15    22    14;
          19    11    20    12;
          17     9    16     8;
          18    10    21    13];

% layerSpecs should be a cell array of ints (1-16) specifying which
% channels in the 16 given by useTheseTrodes to include in each 
% cortical layer
% e.g.,
% layerSpecs{1}=[1 2]; layerSpecs{2}=[3 4 5 6]; etc.
% layerSpecs{1}=[1 2 3];
% layerSpecs{2}=[4 5 6 7];
% layerSpecs{3}=[8 9];
% layerSpecs{4}=[10 11 12 13];
% layerSpecs{5}=[14 15 16];
layerSpecs{1}=1;
layerSpecs{2}=2;
layerSpecs{3}=3;
layerSpecs{4}=4;
layerSpecs{5}=5;
layerSpecs{6}=6;
layerSpecs{7}=7;
layerSpecs{8}=8;
layerSpecs{9}=9;
layerSpecs{10}=10;
layerSpecs{11}=11;
layerSpecs{12}=12;
layerSpecs{13}=13;
layerSpecs{14}=14;
layerSpecs{15}=15;
layerSpecs{16}=16;

% Different LEDconds and time windows for which to calculate trial-by-trial
% FRs across layers
% LEDblocks={3; 3; 3; 3};
% noLEDblocks={1; 1; 1; 1};
% timeWindows={[3.3 4.3];[3.3 3.6];[3.6 3.9];[4 4.3]};
%          spont  ev(minus base) ev
% stimBlocks={[1:2:16]; [2:2:16]; [2:2:16]};
% LEDblocks={[5]; [5]; [5]};
% noLEDblocks={[0]; [0]; [0]};
% timeWindows={[1.25 1.325]; [1.25 1.325]; [1.25 1.325]};
% isevoked={'no'; 'yes'; 'no'};
% spontWindow={[]; [0 1]; []};

% stimBlocks={[1:9]};
% LEDblocks={[5]};
% noLEDblocks={[0]};
% timeWindows={[1.25 1.325]};
% isevoked={'no'};
% spontWindow={[]};

stimBlocks={[1:16];[1:16]};
LEDblocks={[5];[5]};
noLEDblocks={[0];[0]};
timeWindows={[1.325 2.5];[1.325 2.5]};
isevoked={'no';'no'};
spontWindow={[];[]};
stuporBlocks={0; 1};

RigDef=RigDefs();
chSetup=RigDef.ChannelOrder{2}; % This is the 1x16 channel setup
for i=1:length(layerSpecs)
    layerSpecs{i}=chSetup(layerSpecs{i});
end

for i=1:length(layerSpecs)
    disp(i);
    % Calculate PSTHs for each cortical layer individually
    % First get spikes for this cortical layer
    theseChs=layerSpecs{i};
    currSpikes=[];
    for j=1:length(useTheseTrodes)
        if any(ismember(trodeChs(j,:),theseChs))
            evChsDefault=[1 2 3 4];
            includeChs=ismember(trodeChs(j,:),theseChs);
            includeEvChs=evChsDefault(includeChs);
            disp('trode');
            disp(useTheseTrodes(j));
            disp('event_channels');
            disp(includeEvChs);
            spikes=load([RigDef.Dir.Spikes expt.sort.trode(useTheseTrodes(j)).spikesfile]);
            spikes=spikes.spikes;
            spikes.event_channel=spikes.info.detect.event_channel';
            stuporTrials=spikes.sweeps.trials(F1amps'>=0.25);
            noStuporTrials=spikes.sweeps.trials(F1amps'<0.25);
%             spikes.stupor=ismember(spikes.trials,stuporTrials);
            spikes=filtspikes(spikes,0,'fileInd',useFileInds);
%             if length(spikes.stupor)~=length(spikes.led)
%                 disp('error');
%                 return
%             end
%             spikes=filtspikes(spikes,0,'stimcond',useStimConds);
%             spikes=filtspikes(spikes,0,'trigger',2:12);
            spikes=getTheseCh_spikes(spikes,includeEvChs);
            if isempty(currSpikes)
                currSpikes=spikes;
            else
                currSpikes=concatSpikes_noAssigns(currSpikes,spikes);
            end
        end
    end 
    spikes=currSpikes;
    for j=1:length(LEDblocks)
        useStimConds=stimBlocks{j};
        % Now calculate FR across trials using these spikes
        if strcmp(isevoked{j},'no')
            subSpikes=filtspikes(spikes,0,'led',noLEDblocks{j},'stimcond',useStimConds);
            if stuporBlocks{j}
                subSpikes=filtspikes(subSpikes,0,'trials',stuporTrials);
            else
                subSpikes=filtspikes(subSpikes,0,'trials',noStuporTrials);
            end
            [m,s,n]=calcMeanAndStdDuringWindow(subSpikes,timeWindows{j});
            FRs_noLED{j}=n;
            subSpikes=filtspikes(spikes,0,'led',LEDblocks{j},'stimcond',useStimConds);
            if stuporBlocks{j}
                subSpikes=filtspikes(subSpikes,0,'trials',stuporTrials);
            else
                subSpikes=filtspikes(subSpikes,0,'trials',noStuporTrials);
            end
            [m,s,n]=calcMeanAndStdDuringWindow(subSpikes,timeWindows{j});
            FRs_LED{j}=n;
        else
            sWindow=spontWindow{j};
            subSpikes=filtspikes(spikes,0,'led',noLEDblocks{j},'stimcond',useStimConds);
            if stuporBlocks{j}
                subSpikes=filtspikes(subSpikes,0,'trials',stuporTrials);
            else
                subSpikes=filtspikes(subSpikes,0,'trials',noStuporTrials);
            end
            [m,s,n]=calcMeanAndStdDuringWindow(subSpikes,timeWindows{j});
            [spont_m,spont_s,spont_n]=calcMeanAndStdDuringWindow(subSpikes,sWindow);
            FRs_noLED{j}=n-spont_n;
            subSpikes=filtspikes(spikes,0,'led',LEDblocks{j},'stimcond',useStimConds);
            if stuporBlocks{j}
                subSpikes=filtspikes(subSpikes,0,'trials',stuporTrials);
            else
                subSpikes=filtspikes(subSpikes,0,'trials',noStuporTrials);
            end
            [m,s,n]=calcMeanAndStdDuringWindow(subSpikes,timeWindows{j});
            [spont_m,spont_s,spont_n]=calcMeanAndStdDuringWindow(subSpikes,sWindow);
            FRs_LED{j}=n-spont_n;
        end
    end
    layers{i}.FRs_noLED=FRs_noLED;
    layers{i}.FRs_LED=FRs_LED;
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


    
function tryThisToMakePlot
for i=1:16
curr_noLED=layers{i}.FRs_noLED{1};
curr_LED=layers{i}.FRs_LED{1};
noLED_means_spont(i)=mean(curr_noLED);
LED_means_spont(i)=mean(curr_LED);
noLED_stds_spont(i)=std(curr_noLED);
LED_stds_spont(i)=std(curr_LED);
end
figure(); 
plot(1:16,noLED_means_spont,'Color','k');
hold on;
plot(1:16,LED_means_spont,'Color','r');

% for i=1:16
% curr_noLED=layers{i}.FRs_noLED{2}-layers{i}.FRs_noLED{1};
% curr_LED=layers{i}.FRs_LED{2}-layers{i}.FRs_LED{1};
% noLED_means_ev(i)=mean(curr_noLED);
% LED_means_ev(i)=mean(curr_LED);
% noLED_stds_ev(i)=std(curr_noLED);
% LED_stds_ev(i)=std(curr_LED);
% end




% x_1=x1; y_1=y1; cx_1=x; cy_1=y;

% all_control_y=mean([cy_1;cy_2],1); all_expt_y=mean([y_1;y_2],1); all_control_x=cx_1; all_expt_x=x_1;

% [down_control_y,down_expt_y,down_x]=useLargerBinsize(all_control_y,all_expt_y,all_control_x,3);
% figure(); plot(all_control_x,all_control_y,'Color','k'); hold on; plot(all_expt_x,all_expt_y,'Color','r');
% [A,tau,normOfResid]=fitExponential(all_expt_x,all_expt_y,3.3);
% disp(tau);