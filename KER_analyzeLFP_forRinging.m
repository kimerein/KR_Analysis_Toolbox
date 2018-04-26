function [LFPbySweep,Fs,bandPassedLFPbySweep]=KER_analyzeLFP_forRinging(daqFileNames, LFPdata, ONlength, trialLength, chooseNSweeps)
    
    Fs=32000;
    LPcutoff=300;
    HPcutoff=0;
    physChannel=4;
    photoChannel=2;
    DownSampleFactor=10;
    
    if ~isempty(LFPdata)
        LFPbySweep=LFPdata;
        Fs2=Fs/DownSampleFactor;
    else
        LFPbySweep=[];
        for i=1:length(daqFileNames)
            daqFileName=daqFileNames{i};
            [thisLFPbySweep,Fs2,photodiodeData]=KER_readStimTriggeredLFP_withPhotodiode(daqFileName,Fs,LPcutoff,HPcutoff,physChannel,photoChannel);
            LFPbySweep=[LFPbySweep; thisLFPbySweep];
        end
    end
    Fs=Fs2;
    
    sweepInds=randperm(size(LFPbySweep,1));
    LFPbySweep=LFPbySweep(sort(sweepInds(1:chooseNSweeps)),:);
    
    % Align photodiode and LFP response
    trialLength2=size(LFPbySweep,2)*(1/Fs);
    figure;
    %plot(0:trialLength2/(size(LFPbySweep,2)-1):trialLength2,(mean(LFPbySweep,1)-mean(mean(LFPbySweep,1),2))/max(mean(LFPbySweep,1)),'Color','b');
    var1=mean(LFPbySweep,1)/norm(mean(LFPbySweep,1));
    plot(0:trialLength2/(size(LFPbySweep,2)-1):trialLength2,var1-var1(1),'Color','b');
    hold on;
    
    
    %plot(0:trialLength2/(size(LFPbySweep,2)-1):trialLength2,(photodiodeData(1:size(LFPbySweep,2))'-mean(photodiodeData(1:size(LFPbySweep,2))))/max(photodiodeData(1:size(LFPbySweep,2))),'Color','k');
    var2=photodiodeData(1:size(LFPbySweep,2))'/norm(photodiodeData(1:size(LFPbySweep,2))')/100;
    plot(0:trialLength2/(size(LFPbySweep,2)-1):trialLength2,var2-var2(1),'Color','k');
    
   
%     for i=0:size(LFPbySweep,1)-1
%         var2=photodiodeData(i*size(LFPbySweep,2)+1:(i+1)*size(LFPbySweep,2))'/norm(photodiodeData(i*size(LFPbySweep,2)+1:(i+1)*size(LFPbySweep,2))')/100;
%         plot(0:trialLength2/(size(LFPbySweep,2)-1):trialLength2,var2-var2(1),'Color','k');
%     end
    
    
    KER_plotOrganizedLFP(LFPbySweep,Fs,[],1:size(LFPbySweep,1));
    bandPassedLFPbySweep=KER_bandPass_LFP(LFPbySweep,Fs,30,80,1,[],1:size(LFPbySweep,1));
end