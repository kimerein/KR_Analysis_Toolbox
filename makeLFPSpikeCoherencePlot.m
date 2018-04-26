function [allC,allCV1,allC2,xtimes,f,fakeTogetherLGN,fakeTogetherV1]=makeLFPSpikeCoherencePlot(L_dLGN,spikes_dLGN,L_V1,spikes_V1,params,useTrials,setting)

% setting == 1: do LGN spike-LFP coherence plot
allC=[];
allCV1=[];
allC2=[];
xtimes=[];
f=[];
% setting == 2: do LGN and V1 spikes spectrograms
% setting == 3: do LGN and V1 LFP spectrograms
% setting == 4: do LGN and V1 LFP spectrograms
fakeTogetherLGN=[];
fakeTogetherV1=[];

if setting==1
    downSamp=1;
    % movingwin=[2.5 0.01];
    movingwin=[3 0.1];
    trialDuration=10.5;
    freqWindow=[55.5 70];
    
    allC=[];
    for i=1:length(useTrials)
        tr=useTrials(i);
        [C,phi,S12,S1,S2,t,f]=cohgramcpb(downSampAv(L_dLGN(tr,:),downSamp)',downSampAv(spikes_dLGN(tr,1:end-1),downSamp)',movingwin,params);
        allC=[allC; C];
    end
    
    xtimes=linspace(0,length(useTrials)*trialDuration,size(allC,1));
    figure();
    imagesc(xtimes,f(f>=freqWindow(1) & f<=freqWindow(2)),allC(:,f>=freqWindow(1) & f<=freqWindow(2))');
elseif setting==2
%     downSamp=5;
    downSamp=1;
%     movingwin=[2.5 0.01];
%     movingwin=[3 0.1];
    movingwin=[1 0.005];
    trialDuration=10.5;
    freqWindow=[55.5 70];
    params.Fs=params.Fs/downSamp;
%     usewin=[0 0.5];
    usewin=[0 0.5 8 10.5];
%     usewin=[0 10.5];
    x=linspace(0,trialDuration,size(spikes_dLGN,2)-1);
    fakeTogetherLGN=[];
    fakeTogetherV1=[];
    
    allC=[];
    for i=1:length(useTrials)
        tr=useTrials(i);
        curr=spikes_dLGN(tr,1:end-1);
        [C,t,f]=mtspecgrampb(downSampAv([curr(x>=usewin(1) & x<=usewin(2)) curr(x>=usewin(3) & x<=usewin(4))],downSamp)',movingwin,params);
        fakeTogetherLGN=[fakeTogetherLGN curr(x>=usewin(1) & x<=usewin(2)) curr(x>=usewin(3) & x<=usewin(4))];
        allC=[allC; C];
    end
    
    allCV1=[];
    for i=1:length(useTrials)
        tr=useTrials(i);
        curr=spikes_V1(tr,1:end-1);
        [C,t,f]=mtspecgrampb(downSampAv([curr(x>=usewin(1) & x<=usewin(2)) curr(x>=usewin(3) & x<=usewin(4))],downSamp)',movingwin,params);
        fakeTogetherV1=[fakeTogetherV1 curr(x>=usewin(1) & x<=usewin(2)) curr(x>=usewin(3) & x<=usewin(4))];
        allCV1=[allCV1; C];
    end
    
    xtimes=linspace(0,length(useTrials)*((usewin(2)-usewin(1))+(usewin(4)-usewin(3))),size(allC,1));
    figure();
    subplot(2,1,1);
    imagesc(xtimes,f(f>=freqWindow(1) & f<=freqWindow(2)),allC(:,f>=freqWindow(1) & f<=freqWindow(2))');
    subplot(2,1,2);
    imagesc(xtimes,f(f>=freqWindow(1) & f<=freqWindow(2)),allCV1(:,f>=freqWindow(1) & f<=freqWindow(2))');
elseif setting==3
    downSamp=1;
    % movingwin=[2.5 0.01];
%     movingwin=[3 0.1];
    movingwin=[1 0.005];
    trialDuration=10.5;
    freqWindow=[55.5 70];
    params.Fs=params.Fs/downSamp;
    usewin=[0 0.5 8 10.5];
    x=linspace(0,trialDuration,size(L_dLGN,2)-1);
    
    allC=[];
    for i=1:length(useTrials)
        tr=useTrials(i);
        curr=L_dLGN(tr,:);
        [C,t,f]=mtspecgramc(downSampAv([curr(x>=usewin(1) & x<=usewin(2)) curr(x>=usewin(3) & x<=usewin(4))],downSamp)',movingwin,params);
        allC=[allC; C];
    end
    
    allCV1=[];
    for i=1:length(useTrials)
        tr=useTrials(i);
        curr=L_V1(tr,:);
        [C,t,f]=mtspecgramc(downSampAv([curr(x>=usewin(1) & x<=usewin(2)) curr(x>=usewin(3) & x<=usewin(4))],downSamp)',movingwin,params);
        allCV1=[allCV1; C];
    end
    
%     xtimes=linspace(0,length(useTrials)*trialDuration,size(allC,1));
    xtimes=linspace(0,length(useTrials)*((usewin(2)-usewin(1))+(usewin(4)-usewin(3))),size(allC,1));
    figure();
    subplot(2,1,1);
    imagesc(xtimes,f(f>=freqWindow(1) & f<=freqWindow(2)),allC(:,f>=freqWindow(1) & f<=freqWindow(2))');
    subplot(2,1,2);
    imagesc(xtimes,f(f>=freqWindow(1) & f<=freqWindow(2)),allCV1(:,f>=freqWindow(1) & f<=freqWindow(2))');
elseif setting==4
    downSamp=5;
    % movingwin=[2.5 0.01];
    movingwin=[10 0.5];
    trialDuration=10.5;
    freqWindow=[55.5 70];
    params.Fs=params.Fs/downSamp;
    
    allC=[];
    for i=1:length(useTrials)
        tr=useTrials(i);
        temp1=downSampMatrix(L_dLGN,downSamp);
        temp2=downSampMatrix(spikes_V1(:,1:end-1),downSamp);
        [C,phi,S12,S1,S2,t,f]=cohgramcpb(temp1(tr,:)',temp2(tr,:)',movingwin,params);
        allC=[allC; C];
    end
    
    allCV1=[];
    for i=1:length(useTrials)
        tr=useTrials(i);
        temp1=downSampMatrix(L_V1,downSamp);
        temp2=downSampMatrix(spikes_dLGN(:,1:end-1),downSamp);
        [C,phi,S12,S1,S2,t,f]=cohgramcpb(temp1(tr,:)',temp2(tr,:)',movingwin,params);
        allCV1=[allCV1; C];
    end
  
    allC2=[];
    for i=1:length(useTrials)
        tr=useTrials(i);
        temp1=downSampMatrix(spikes_dLGN(:,1:end-1),downSamp);
        temp2=downSampMatrix(spikes_V1(:,1:end-1),downSamp);
        [C,phi,S12,S1,S2,t,f]=cohgrampb(temp1(tr,:)',temp2(tr,:)',movingwin,params);
        allC2=[allC2; C];
    end
    
    xtimes=linspace(0,length(useTrials)*trialDuration,size(allC,1));
    figure();
    subplot(3,1,1);
    imagesc(xtimes,f(f>=freqWindow(1) & f<=freqWindow(2)),allC(:,f>=freqWindow(1) & f<=freqWindow(2))');
    subplot(3,1,2);
    imagesc(xtimes,f(f>=freqWindow(1) & f<=freqWindow(2)),allCV1(:,f>=freqWindow(1) & f<=freqWindow(2))');
    subplot(3,1,3);
    imagesc(xtimes,f(f>=freqWindow(1) & f<=freqWindow(2)),allC2(:,f>=freqWindow(1) & f<=freqWindow(2))');
end
    