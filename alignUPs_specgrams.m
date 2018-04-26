function [finalAlignedSpecgram,specgrams]=alignUPs_specgrams(UPstates,trialLEDs,useLED,trialStims,useStimcond,LEDwindow,specgram_x,specgrams)

trialDuration=5;
dodownSamp=1;
downSampFactor=10;

if isempty(specgram_x)
    a=size(specgrams{1});
    specgram_x=linspace(0,trialDuration,a(2));
end
if dodownSamp==1
    specgram_x=downSampAv(specgram_x,downSampFactor);
    for i=1:length(specgrams)
        currSpec=specgrams{i};
        disp(i);
        currSpec=currSpec';
        for k=1:size(currSpec,2)
            currSpec(:,k)=smooth(currSpec(:,k),downSampFactor);
        end
        currSpec=downsample(currSpec,downSampFactor);
        specgrams{i}=currSpec';
    end
end

currLine=cell(1,length(specgrams));
for j=1:a(1) % Iterate through lines of specgram
    disp(j);
    for i=1:length(specgrams)
        currSpec=specgrams{i};
        currLine{i}=currSpec(j,:);
    end
    [~,~,xled,yled]=alignUPs_general(UPstates,trialLEDs,useLED,trialStims,useStimcond,LEDwindow,specgram_x,currLine);
    finalAlignedSpecgram(j,:)=yled;
end

figure(); 
imagesc(finalAlignedSpecgram);

    
        