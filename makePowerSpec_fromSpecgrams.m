function [currTrialSpecgrams,x,av,lowBound,highBound,LFPspecgram]=makePowerSpec_fromSpecgrams(LFPbySweep,Fs,alreadyTrialSpecgrams)

dataDir='Z:\Updating Data\Analysis Computer\All Specgram Data\Mawake54\Thin Time Bins\';
fileStepSize=5;
startFile=94;
endFile=133;
theseFiles=startFile:fileStepSize+1:endFile;

if ~isempty(alreadyTrialSpecgrams)
    currTrialSpecgrams=alreadyTrialSpecgrams;
else
    currTrialSpecgrams=[];
    for i=1:length(theseFiles)-1
        currFiles=[theseFiles(i) theseFiles(i+1)-1];
        for k=1:1+fileStepSize
            currDetails=['fileAIs_' num2str(currFiles(1)) 'to' num2str(currFiles(2)) '_details_' num2str(k) '.mat'];
            currSpecgrams=['fileAIs_' num2str(currFiles(1)) 'to' num2str(currFiles(2)) '_specgrams_' num2str(k) '.mat'];
            disp(['Loading' currSpecgrams]);
            temp=load([dataDir currDetails]);
            trialDetails=temp.trialDetails;
            trialSpecgram.otherDetails=trialDetails.otherDetails;
            trialSpecgram.ledConds=trialDetails.ledConds;
            trialSpecgram.stimConds=trialDetails.stimConds;
            temp=load([dataDir currSpecgrams]);
            trialSpecgram.specgrams=temp.a;
            clear temp;
            currTrialSpecgrams=bootstrapPowerSpec_fromSpecgrams(currTrialSpecgrams,trialSpecgram,LFPbySweep,Fs,0);
        end
    end
end
    
[~,x,av,lowBound,highBound,LFPspecgram]=bootstrapPowerSpec_fromSpecgrams([],currTrialSpecgrams,LFPbySweep,Fs,1);
        
        
        