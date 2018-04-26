function tryDirectConvolution(name,spikes,xtimes,timeCourse,transferFunction,dp)

shiftBy3ms=1;

saveDir='W:\Analysis Computer\Across Mice Finals\Frequency Response\';
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
Fs=32000;
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
params.useFileInd=[1:10000];
params.binsize=1; % in ms - MUST MATCH xtimes units
params.useStimcond=[1:128];
params.useLEDcond=[1.050 2.050 4.050 6.050 8.050 10.050 12.050 14.050 16.050 18.050 20.050 30.050 40.050 50.050 60.050];
params.trialDuration=4; % in s
params.stimulusOn=[1 3]; % in s relative to trial onset
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
if isempty(dp)
    save([saveDirName '\params.mat'],'params');

    spikes=filtspikes(spikes,0,'fileInd',params.useFileInd);
    criterionSpikes=filtspikes(spikes,0,'stimcond',params.useStimcond);
    ma=max(timeCourse);
    scaleFactor=1/ma;
    shutoff=timeCourse.*scaleFactor;
    for i=1:length(params.useLEDcond)
        disp('making predictions');
        ledValue=params.useLEDcond(i);
        temp=[];
        temp1=[];
        for j=1:length(ledValue)
            criterionSpikes=makeTempField(criterionSpikes,'led',ledValue(j));
            temp(j,:)=criterionSpikes.temp;
            temp1(j,:)=criterionSpikes.sweeps.temp;
            criterionSpikes.temp=[];
            criterionSpikes.sweeps.temp=[];
        end
        criterionSpikes.temp=sum(temp,1)>=1;
        criterionSpikes.sweeps.temp=sum(temp1,1)>=1;
        currSpikes=filtspikes(criterionSpikes,0,'temp',1);
        [~,~,~,xpoints1,ypoints1]=psth_wStdev_valuesOnly(currSpikes,params.binsize);
        if shiftBy3ms==1
%             output=conv(shutoff,ypoints1); 
            output=conv(ypoints1,shutoff);
            output=[zeros(1,3) output];
        else
            output=conv(ypoints1,shutoff);     
        end
%         [result,x]=getConvForFreq(ypoints1,0,xtimes,timeCourse,transferFunction);
        directPrediction.x{i}=xpoints1;
        directPrediction.input{i}=ypoints1;
        directPrediction.result{i}=output(1:length(xpoints1));
    end
    save([saveDirName '\directPrediction.mat'],'directPrediction');
else
    directPrediction=dp;
end
coh=zeros(1,length(params.useLEDcond));
% cerr1=zeros(1,length(params.useLEDcond));
% cerr2=zeros(1,length(params.useLEDcond));
for i=1:length(params.useLEDcond)
    ledValue=params.useLEDcond(i);
    currdp=directPrediction.result{i};
    currx=directPrediction.x{i};
    stimWindow=[1.3 3];
    [c]=calcCoherenceContinuous(freqs(i),currdp(currx>=stimWindow(1) & currx<=stimWindow(2)),stimWindow);
    coh(i)=c;
%     cerr1(i)=c1;
%     cerr2(i)=c2;
end
cohereDP.coh=coh;
cohereDP.freqs=freqs;
% cohereDP.cerr1=cerr1;
% cohereDP.cerr2=cerr2;
save([saveDirName '\cohereDP.mat'],'cohereDP');
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