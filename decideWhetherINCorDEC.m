function [INCresponses,DECresponses,INCFarresponses,DECFarresponses,x,ampAnalysis,phases,popTrialData]=decideWhetherINCorDEC(spikes,useAssigns,name)

skip=1;

useS=[1:8];
% useS=1:32;
bin=1;
nStdevs=1;
useLED=0;
baseWindow=[0 1];
stimWindow=[1 4];
analyzeWindow=[4 5];
params.tapers=[10 18];
params.pad=0;
params.fpass=[1 70];
params.Fs=1000/bin;
params.trialave=1;   

phases=[];
if skip==0
    [x,psths,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED]=getTrialByTrialUnitPSTH(spikes,useAssigns);
    x=[x x(end)+(x(2)-x(1))];
    p=psths{1};
    INCresponses=[];
    DECresponses=[];
    INCFarresponses=[];
    DECFarresponses=[];
    unitClassification=nan(length(psths),length(useS));
    for i=1:length(psths)
        p=psths{i};
        allBase=nanmean(nanmean(p(:,x>=baseWindow(1) & x<=baseWindow(2)),1),2);
        baseStd=nanstd(downSampAv(nanmean(p(:,x>=baseWindow(1) & x<=baseWindow(2)),1),50));
        for j=1:length(useS)
            currS=useS(j);
            stimResponse=nanmean(p(ismember(unitByUnitLED{i},useLED) & ismember(unitByUnitStimcond{i},currS),:),1);
            stimRAv=nanmean(stimResponse(x>=stimWindow(1) & x<=stimWindow(2)));
            if stimRAv>allBase
                INCresponses=[INCresponses; stimResponse];
                unitClassification(i,j)=1;
            elseif stimRAv<allBase
                DECresponses=[DECresponses; stimResponse];
                unitClassification(i,j)=2;
            end
            if stimRAv>=nStdevs*baseStd+allBase
                INCFarresponses=[INCFarresponses; stimResponse];
                unitClassification(i,j)=1.5;
            elseif stimRAv<=allBase-nStdevs*baseStd
                DECFarresponses=[DECFarresponses; stimResponse];
                unitClassification(i,j)=2.5;
            end
        end
    end
    
    trialINC=zeros(size(psths{1}));
    trialDEC=zeros(size(psths{1}));
    trialINCFar=zeros(size(psths{1}));
    trialDECFar=zeros(size(psths{1}));
    trialNOMOD=zeros(size(psths{1}));
    % Get trial by trial INCs and DECs
    for i=1:length(psths)
        p=psths{i};
        for j=1:length(useS)
            currS=useS(j);
            useTrials=ismember(unitByUnitLED{i},useLED) & ismember(unitByUnitStimcond{i},currS);
            if unitClassification(i,j)==1 || unitClassification(i,j)==1.5
                trialINC(useTrials,:)=trialINC(useTrials,:)+p(useTrials,:);
            elseif unitClassification(i,j)==2 || unitClassification(i,j)==2.5
                trialDEC(useTrials,:)=trialDEC(useTrials,:)+p(useTrials,:);
            end
            if unitClassification(i,j)==1.5
                trialINCFar(useTrials,:)=trialINCFar(useTrials,:)+p(useTrials,:);
            elseif unitClassification(i,j)==2.5
                trialDECFar(useTrials,:)=trialDECFar(useTrials,:)+p(useTrials,:);
            end
            if unitClassification(i,j)==1 || unitClassification(i,j)==2
                trialNOMOD(useTrials,:)=trialNOMOD(useTrials,:)+p(useTrials,:);
            end
        end
    end
    
    % Find phase and amplitude over analyzeWindow
    % Use sliding window
    % for trialINC
    [S,t,f]=mtspecgrampb(trialINC',[1.5 0.05],params);
    [S_DEC,t_DEC,f_DEC]=mtspecgrampb(trialDEC',[1.5 0.05],params);
    [S_INCFar,t_INCFar,f_INCFar]=mtspecgrampb(trialINCFar',[1.5 0.05],params);
    [S_DECFar,t_DECFar,f_DECFar]=mtspecgrampb(trialDECFar',[1.5 0.05],params);
    [S_NOMOD,t_NOMOD,f_NOMOD]=mtspecgrampb(trialNOMOD',[1.5 0.05],params);
    ampAnalysis.S=S;
    ampAnalysis.t=t;
    ampAnalysis.f=f;
    ampAnalysis.S_DEC=S_DEC;
    ampAnalysis.S_INCFar=S_INCFar;
    ampAnalysis.S_DECFar=S_DECFar;
    ampAnalysis.S_NOMOD=S_NOMOD;
    save(['W:\Analysis Computer\150421 Alpha\' name 'nothing.mat']);
else
    load('W:\Analysis Computer\150421 Alpha\stateOfProgram1_backup.mat');
end

trialNOMOD=zeros(size(psths{1}));
for i=1:length(psths)
    p=psths{i};
    for j=1:length(useS)
        currS=useS(j);
        useTrials=ismember(unitByUnitLED{i},useLED) & ismember(unitByUnitStimcond{i},currS);
        if unitClassification(i,j)==1 || unitClassification(i,j)==2
            trialNOMOD(useTrials,:)=trialNOMOD(useTrials,:)+p(useTrials,:);
        end
    end
end
[S_NOMOD,t_NOMOD,f_NOMOD]=mtspecgrampb(trialNOMOD',[1.5 0.05],params);
ampAnalysis.S_NOMOD=S_NOMOD;
popTrialData.trialNOMOD=trialNOMOD;
return
popTrialData.trialINC=trialINC;
popTrialData.trialDEC=trialDEC;
popTrialData.trialINCFar=trialINCFar;
popTrialData.trialDECFar=trialDECFar;
% find frequency with max amplitude between 4 and 8 Hz at offset
S_analyzeWindow=nanmean([S(t>=analyzeWindow(1) & t<=analyzeWindow(2),:); S_DEC(t>=analyzeWindow(1) & t<=analyzeWindow(2),:); S_INCFar(t>=analyzeWindow(1) & t<=analyzeWindow(2),:); S_DECFar(t>=analyzeWindow(1) & t<=analyzeWindow(2),:)],1);
subS=S_analyzeWindow(f>=4 & f<=10);
subf=f(f>=4 & f<=10);
[~,mind]=max(subS);
maxAmpFreq=subf(mind);
% make that frequency for looking at phase of spikes
fakex=0:bin/1000:(analyzeWindow(2)-analyzeWindow(1))-bin/1000;
fakey=sin(2*pi*maxAmpFreq.*fakex);
% get phase relative to this frequency at offset for each population
[C,phi,~,~,~,f_coh]=coherencypb(repmat(fakey',1,size(trialINC,1)),trialINC(:,x>=analyzeWindow(1) & x<=analyzeWindow(2))',params,0);
phases.INC=nanmean(phi([find(f_coh>=maxAmpFreq,1,'first') find(f_coh<=maxAmpFreq,1,'last')]));
[C,phi,~,~,~,f_coh]=coherencypb(repmat(fakey',1,size(trialINC,1)),trialDEC(:,x>=analyzeWindow(1) & x<=analyzeWindow(2))',params,0);
phases.DEC=nanmean(phi([find(f_coh>=maxAmpFreq,1,'first') find(f_coh<=maxAmpFreq,1,'last')]));
[C,phi,~,~,~,f_coh]=coherencypb(repmat(fakey',1,size(trialINC,1)),trialINCFar(:,x>=analyzeWindow(1) & x<=analyzeWindow(2))',params,0);
phases.INCFar=nanmean(phi([find(f_coh>=maxAmpFreq,1,'first') find(f_coh<=maxAmpFreq,1,'last')]));
[C,phi,~,~,~,f_coh]=coherencypb(repmat(fakey',1,size(trialINC,1)),trialDECFar(:,x>=analyzeWindow(1) & x<=analyzeWindow(2))',params,0);
phases.DECFar=nanmean(phi([find(f_coh>=maxAmpFreq,1,'first') find(f_coh<=maxAmpFreq,1,'last')]));

% Get phases of populations relative to each other
[C,phi,~,~,~,f_coh]=coherencypb(trialINC(:,x>=analyzeWindow(1) & x<=analyzeWindow(2))',trialDEC(:,x>=analyzeWindow(1) & x<=analyzeWindow(2))',params,0);
phases.INCvsDEC=nanmean(phi([find(f_coh>=maxAmpFreq,1,'first') find(f_coh<=maxAmpFreq,1,'last')]));
[C,phi,~,~,~,f_coh]=coherencypb(trialINCFar(:,x>=analyzeWindow(1) & x<=analyzeWindow(2))',trialDECFar(:,x>=analyzeWindow(1) & x<=analyzeWindow(2))',params,0);
phases.INCFarvsDECFar=nanmean(phi([find(f_coh>=maxAmpFreq,1,'first') find(f_coh<=maxAmpFreq,1,'last')]));


        
        
        