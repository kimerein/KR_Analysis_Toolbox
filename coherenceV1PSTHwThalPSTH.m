function coherenceV1PSTHwThalPSTH(name,V1spikes,thalspikes,psths)

% saveDir='W:\Analysis Computer\Across Mice Finals\Frequency Response\';
saveDir='W:\Analysis Computer\Across Mice Finals\Coherence\';
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
Fs=32000;
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
params.useFileInd=[];
params.binsize=1; % in ms - MUST MATCH xtimes units
params.useStimcond=[1:128];
params.useLEDcond=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% params.useLEDcond=[1.050 2.050 4.050 6.050 8.050 10.050 12.050 14.050 16.050 18.050 20.050 30.050 40.050 50.050 60.050];
params.trialDuration=4; % in s
stimWindow=[1.3 3];
params.stimulusOn=stimWindow; % in s relative to trial onset
params.stimulusType='Sine';
params.LEDonset=500; % ms into trial for LED onset
params.LEDduration=3000; % in ms, duration of LED pulse -- maybe short to exclude photo-artifact
params.baseWindow=[0 1];
params.analysisType='cross-corr'; 
params.anesthOrAwake='anesth';
params.anesthType='iso';
params.LEDintensity=5;
% Set up save directory for this mouse
saveDirName=[saveDir name];
if exist(saveDirName,'dir')==0
    st=mkdir(saveDirName);
    if ~st
        disp('Could not create save directory');
    end
end

save([saveDirName '\params.mat'],'params');

% V1spikes=filtspikes(V1spikes,0,'fileInd',params.useFileInd);
% thalspikes=filtspikes(thalspikes,0,'fileInd',params.useFileInd);
if isempty(psths)
V1spikes=filtspikes(V1spikes,0,'stimcond',params.useStimcond);
thalspikes=filtspikes(thalspikes,0,'stimcond',params.useStimcond);
for i=1:length(params.useLEDcond)
    disp('calculating PSTHs');
    ledValue=params.useLEDcond(i);
    temp=[];
    temp1=[];
    for j=1:length(ledValue)
        V1spikes=makeTempField(V1spikes,'led',ledValue(j));
        temp(j,:)=V1spikes.temp;
        temp1(j,:)=V1spikes.sweeps.temp;
        V1spikes.temp=[];
        V1spikes.sweeps.temp=[];
    end
    V1spikes.temp=sum(temp,1)>=1;
    V1spikes.sweeps.temp=sum(temp1,1)>=1;
    currV1Spikes=filtspikes(V1spikes,0,'temp',1);
    temp=[];
    temp1=[];
    for j=1:length(ledValue)
        thalspikes=makeTempField(thalspikes,'led',ledValue(j));
        temp(j,:)=thalspikes.temp;
        temp1(j,:)=thalspikes.sweeps.temp;
        thalspikes.temp=[];
        thalspikes.sweeps.temp=[];
    end
    thalspikes.temp=sum(temp,1)>=1;
    thalspikes.sweeps.temp=sum(temp1,1)>=1;
    currThalSpikes=filtspikes(thalspikes,0,'temp',1);
    [~,~,~,xpoints1,ypoints1]=psth_wStdev_valuesOnly(currV1Spikes,params.binsize);
    [~,~,~,xpoints1thal,ypoints1thal]=psth_wStdev_valuesOnly(currThalSpikes,params.binsize);
%     output=conv(shutoff,ypoints1);
    %         [result,x]=getConvForFreq(ypoints1,0,xtimes,timeCourse,transferFunction);
    psths.V1_x{i}=xpoints1;
    psths.V1_y{i}=ypoints1;
    psths.thal_x{i}=xpoints1thal;
    psths.thal_y{i}=ypoints1thal;
end
end
save([saveDirName '\psths.mat'],'psths');

coh=zeros(1,length(params.useLEDcond));
phase=zeros(1,length(params.useLEDcond));
% cerr1=zeros(1,length(params.useLEDcond));
% cerr2=zeros(1,length(params.useLEDcond));
for i=1:length(params.useLEDcond)
    ledValue=params.useLEDcond(i);
    currv1_y=psths.V1_y{i};
    currv1_x=psths.V1_x{i};
    currthal_y=psths.thal_y{i};
    currthal_x=psths.thal_x{i};
%     [c,p]=calcCoherenceTwoBinnedPP(freqs(i),currv1_y(currv1_x>=stimWindow(1) & currv1_x<=stimWindow(2)),currthal_y(currthal_x>=stimWindow(1) & currthal_x<=stimWindow(2)));
    [c,p]=calcCoherenceTwoBinnedPP(freqs(i),currthal_y(currthal_x>=stimWindow(1) & currthal_x<=stimWindow(2)),currv1_y(currv1_x>=stimWindow(1) & currv1_x<=stimWindow(2)));
    coh(i)=c;
    phase(i)=p;
%     cerr1(i)=c1;
%     cerr2(i)=c2;
end
coherePSTH.coh=coh;
coherePSTH.freqs=freqs;
coherePSTH.phase=phase;
% cohereDP.cerr1=cerr1;
% cohereDP.cerr2=cerr2;
save([saveDirName '\coherePSTH.mat'],'coherePSTH');
end

function [out,x]=applyTransferFunction(x,in,transferFunction)
    th=transferFunction.thalamus;
    cx=transferFunction.cortex;
    out=zeros(size(in));
    for i=1:length(in)
        if in(i)<min(th)
            out(i)=0;
        elseif in(i)>max(th)
            out(i)=max(cx);
        else
            if isempty(find(th<in(i),1,'last'))
                out(i)=min(cx);
            elseif isempty(find(th>in(i),1,'first'))
                out(i)=max(cx);
            else
                out(i)=mean([cx(find(th<in(i),1,'last')) cx(find(th>in(i),1,'first'))]);
                if isnan(out(i))
                    disp('help');
                end
            end
        end
    end
end

function [output,x]=getConvForFreq(input,doPlot,t,timeCourse,transferFunction)
if isempty(timeCourse)
    x=0:0.0001:2; % in s
%     x=0:0.0001:2.0031; % in s
%     net_tau=0.012; % in s
    net_tau=0.02; % in s
    shutoff=[ones(1,length(x(x>=0 & x<=0.003))) exp(-x./net_tau)];
%     shutoff=exp(-x./net_tau);
else
    if max(t)<4
        tsteps=t(2)-t(1);
        x=[t max(t)+tsteps:tsteps:4];
        timeCourse=[timeCourse zeros(1,length(max(t)+tsteps:tsteps:4))];
    else
        x=t; % in s
    end
    % Scale starting value to 1
    ma=max(timeCourse);
    scaleFactor=1/ma;
    shutoff=timeCourse.*scaleFactor;
end

if doPlot==1
    figure(); 
    if isempty(timeCourse)
        plot([x 2.0001:0.0001:2+0.003+0.0001],shutoff);
    else
        plot(x,shutoff);
    end
end

% Apply thalamic filter
if doPlot==1
    figure();
    plot(x,input);
    title('Thalamic Response');
end
    
if ~isempty(transferFunction)
    [input,x]=applyTransferFunction(x,input,transferFunction);
end

if doPlot==1
    figure();
    plot(x,input);
    title('After Applying Transfer Function');
end
    
% output=conv(shutoff,input,'same');
output=conv(shutoff,input);
if isempty(timeCourse)
    x=[x 2.0001:0.0001:2+0.003+0.0001];
end

% figure();
% plot(0:x(2)-x(1):(x(2)-x(1))*(length(output)-1),output);
% title('After Convolution with V1 Tau');

end