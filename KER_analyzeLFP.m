function [LFPbySweep,Fs,bandPassedLFPbySweep]=KER_analyzeLFP(spikes, daqFileNames, LFPdata, redTrials, blueTrials, ONlength, ONstart, trialLength)
    
    Fs=32000;
    LPcutoff=300;
    HPcutoff=0;
    physChannel=4;
    DownSampleFactor=10;
    
%     redTrials=[];
%     blueTrials=[];
%     countTrialByTrialOFFspikes=zeros(max(spikes.trials),1);
%     for i=1:max(spikes.trials)
%         someSpikes=filtspikes(spikes,0,'trials',i);
%         countTrialByTrialOFFspikes(i)=length(someSpikes.spiketimes(someSpikes.spiketimes>ONstart+ONlength & someSpikes.spiketimes<=trialLength));
%         if countTrialByTrialOFFspikes(i)>=redThresh
%             %disp(i);
%             redTrials=[redTrials i];
%         else
%             blueTrials=[blueTrials i];
%         end
%     end
%     [h,x]=hist(countTrialByTrialOFFspikes,8);
%     figure;
%     bar(x,h);
    
    KER_raster(spikes,redTrials);
    if ~isempty(LFPdata)
        LFPbySweep=LFPdata;
        Fs2=Fs/DownSampleFactor;
    else
        LFPbySweep=[];
        for i=1:length(daqFileNames)
            daqFileName=daqFileNames{i};
            [thisLFPbySweep,Fs2]=KER_readStimTriggeredLFP(daqFileName,Fs,LPcutoff,HPcutoff,physChannel);
            LFPbySweep=[LFPbySweep; thisLFPbySweep];
        end
    end
    Fs=Fs2;
    
    KER_plotOrganizedLFP(LFPbySweep,Fs,redTrials,blueTrials);
    bandPassedLFPbySweep=KER_bandPass_LFP(LFPbySweep,Fs,30,80,1,redTrials,blueTrials);
    
end