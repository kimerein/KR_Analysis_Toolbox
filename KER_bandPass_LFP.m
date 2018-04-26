function bandPassedLFP=KER_bandPass_LFP(LFPdata,Fs,lowCutoff,highCutoff,compare,redTrials,blueTrials)
% Band-pass filters LFPdata between lowCutoff and highCutoff; plots results
% Compares redTrials and blueTrials, looks for significant differences in
% LFP
% 
% LFPdata is LFP data organized by sweeps (sweeps are rows, different
% samples are columns)
% 
% Fs is current sampling frequency of LFPdata
% 
% lowCutoff is lower cut-off for band-pass filter in Hz
% highCutoff is upper cut-off for band-pass filter in Hz
%
% compare is whether or not to compare trialSet1 with trialSet2
% compare==1: red (redTrials) and blue (blueTrials) LFP traces are plotted
% Locations of significant difference are indicated by black line
% 
% redTrials are the trial numbers in the first set
% blueTrials are the trial numbers in the second set

trialLength=size(LFPdata,2)*(1/Fs);

% Low-pass-filter data
disp('LP-filtering');
LFPdata=fftFilter(LFPdata',Fs,highCutoff,1);
LFPdata=LFPdata';
disp('Done LP-filtering');

% High-pass-filter data
disp('HP-filtering');
LFPdata=fftFilter(LFPdata',Fs,lowCutoff,2);
LFPdata=LFPdata';
disp('Done HP-filtering');

% Plot band-pass-filtered data
% Align all LFP traces to initial value = 0
for i=1:size(LFPdata,1)
    initVal=LFPdata(i,1);
    LFPdata(i,:)=LFPdata(i,:)-initVal;
end

figure;
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata,1),'Color','k');
title('Average Stimulus-Triggered LFP - Initial Values Aligned');

bandPassedLFP=LFPdata;

if compare==1
    figure;
    redTrials=redTrials';
    blueTrials=blueTrials';
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:),1),'Color','r');
    hold on;
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:),1),'Color','b');
    title('Average Stimulus-Triggered LFP Traces for Different Trial Sets - Initial Values Aligned');
    
    redMean=mean(mean(LFPdata(redTrials,:),1),2);
    line([0 trialLength],[redMean redMean],'Color','r');
    blueMean=mean(mean(LFPdata(blueTrials,:),1),2);
    line([0 trialLength],[blueMean blueMean],'Color','b');
    redGammaAmplitude=mean(sum(abs(LFPdata(redTrials,:)-redMean)),1);
    blueGammaAmplitude=mean(sum(abs(LFPdata(blueTrials,:)-blueMean)),1);
    disp(redGammaAmplitude);
    disp(blueGammaAmplitude);
    
    redMeans=zeros(size(LFPdata(redTrials,:),1),size(LFPdata(:,:),2));
    blueMeans=zeros(size(LFPdata(blueTrials,:),1),size(LFPdata(:,:),2));
    rM=mean(LFPdata(redTrials,:),2);
    bM=mean(LFPdata(blueTrials,:),2);
    for i=1:size(LFPdata(redTrials,:),2)
        redMeans(:,i)=rM;
        blueMeans(:,i)=bM;
    end
    
    figure;
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:)-redMeans,1),'Color','r');
    hold on;
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:)-blueMeans,1),'Color','b');
    title('Average Stimulus-Triggered LFP Traces with Mean and Stdev');
    
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:)-redMeans,1)+std(LFPdata(redTrials,:)-redMeans,1),'Color','m');
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:)-redMeans,1)-std(LFPdata(redTrials,:)-redMeans,1),'Color','m');
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:)-blueMeans,1)+std(LFPdata(blueTrials,:)-blueMeans,1),'Color','c');
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:)-blueMeans,1)-std(LFPdata(blueTrials,:)-blueMeans,1),'Color','c');
    statLine=min(mean(LFPdata(redTrials,:)-redMeans,1));
    
    figure;
%     plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:)-redMeans,1),'Color','r');
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:)-redMeans,1),'Color','r');
    hold on;
%     plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:)-blueMeans,1),'Color','b');
    plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:)-blueMeans,1),'Color','b');
    windowSize=0.1;
    for i=windowSize:windowSize:trialLength
        colInds=floor((i-windowSize)/(trialLength/(size(LFPdata,2)-1)))+1:floor(i/(trialLength/(size(LFPdata,2)-1)))+1;
        redMeans=zeros(size(LFPdata(redTrials,colInds),1),length(colInds));
        blueMeans=zeros(size(LFPdata(blueTrials,colInds),1),length(colInds));
        rM=mean(LFPdata(redTrials,colInds),2);
        bM=mean(LFPdata(blueTrials,colInds),2);
        for j=1:length(colInds)
            redMeans(:,j)=rM;
            blueMeans(:,j)=bM;
        end
        segLFPDataRed=sum(abs(LFPdata(redTrials,colInds)-redMeans),2);
        segLFPDataBlue=sum(abs(LFPdata(blueTrials,colInds)-blueMeans),2);
        [h,p,ci]=ttest2(segLFPDataRed,segLFPDataBlue);
        if h==1
             disp(h);
             disp(p);
             disp(ci);
            if mean(segLFPDataRed)>mean(segLFPDataBlue)
                line([i-windowSize; i],[statLine; statLine],'Color','r');
            else
                line([i-windowSize; i],[statLine; statLine],'Color','b');
            end
        end
    end
    disp('trial length = ');
    disp((1/Fs)*size(LFPdata(:,:),2));
end
