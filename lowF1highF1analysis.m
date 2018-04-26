function lowF1highF1analysis(outputDir,dLGNpsth,usel_noLED,usel_LED,noThetaTrials,uses,uses_tri)

% Mawake322
% uses=   {{1 1 1}; {1 1}};
% uses_tri={{[13:72]-12 [141:212]-12 [277:316]-12}; {[73:140]-12 [213:276]-12}};
% Mawake321
% uses=   {{1 1}; {1}};
% uses_tri={{[1:56] [157:200]}; {[57:156]}};
% Mawake320
% uses=   {{1 1}; {1 1}};
% uses_tri={{[1:80] [157:260]}; {[81:156] [261:340]}};
% Mawake319
% uses=   {{1 2 1}; {2 2}; {1}};
% uses_tri={{[1:80] [81:204] [273:332]}; {[1:80] [273:332]}; {[81:204]}};
% Mawake318
% uses=   {{1 1 1}; {2 2}; {2}};
% uses_tri={{[1:80] [81:196] [197:276]}; {[1:80] [197:276]}; {[81:196]}};
% Mawake317
% uses=   {{1 1}; {2 3}; {2}};
% uses_tri={{[1:188] [189:316 349:392]}; {[1:188] [189:316 349:392]}; {[189:316 349:392]}};
% Mawake316
% uses=   {{1 1}; {2}; {2}};
% uses_tri={{[1:136] [189:268]}; {[1:136]}; {[189:268]}};
% Mawake315
% uses=   {{1:2 1}; {3}; {2}};
% uses_tri={{[9:164] [165:224]}; {[9:164]}; {[165:224]}};
% Mawake314
% uses=   {{1:2 1 1}; {3}; {2}; {2}};
% uses_tri={{[1:160] [197:268] [269:308]}; {[1:160]}; {[197:268]}; {[269:308]}};
% Mawake313
% uses=   {{1:2 1}; {3}; {2}};
% uses_tri={{[1:128] [169:244]}; {[1:128]}; {[169:244]}};
% Mawake312
% uses=   {{1:2 1}; {3}; {2}};
% uses_tri={{[33:132]-32 [169:232]-32}; {[33:132]-32}; {[169:232]-32}};
% Mawake311
% uses=   {{1:2 1 1}; {3}; {2}; {2}};
% uses_tri={{[13:148]-12 [204:254]-12 [255:299]-12}; {[13:148]-12}; {[204:254]-12}; {[255:299]-12}};
% Mawake310
% uses=   {{1:2 1:2 1}; {3 3}; {2}};
% uses_tri={{[9:144]-8 [201:228]-8 [229:279]-8}; {[9:144]-8 [201:228]-8}; {[229:279]-8}};
% Mawake206
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}; {10}; {11}; {12}};
% uses_tri={{[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}};

% Mawake326
% uses=   {{1:2}};
% uses_tri={{[1:10000]}};


%uses=unique(dLGNpsth.unitStimcond{1});
usealls=1:10000;

% if isempty(LFPdata)
%     if isempty(LFP_Fs)
%         LFP_Fs=2500;
%     end
%     
%     disp('1');
%     [~,~,LFPbySweep,ledConds,stimConds]=makeCSD2(expt,useFileInd,[],[],[]);
%     save([outputDir '\' 'LFPbySweep.mat'],'LFPbySweep'); % V1 then hippo
%     save([outputDir '\' 'LFPbySweep_ledConds.mat'],'ledConds');
%     save([outputDir '\' 'LFPbySweep_stimConds.mat'],'stimConds');
%     
%     L=LFPbySweep{1};
%     L=L(~isnan(ledConds) & ~isnan(stimConds),:);
%     LFPbySweep{1}=L;
%     L=LFPbySweep{2};
%     L=L(~isnan(ledConds) & ~isnan(stimConds),:);
%     LFPbySweep{2}=L;
% else
%     LFPbySweep{1}=LFPdata.LFPbySweep{1};
%     if length(LFPdata.LFPbySweep)>1
%         LFPbySweep{2}=LFPdata.LFPbySweep{2};
%     end
%     ledConds=LFPdata.ledConds;
%     stimConds=LFPdata.stimConds;
%     
%     L=LFPbySweep{1};
%     L=L(~isnan(ledConds) & ~isnan(stimConds),:);
%     LFPbySweep{1}=L;
%     if length(LFPdata.LFPbySweep)>1
%         L=LFPbySweep{2};
%         L=L(~isnan(ledConds) & ~isnan(stimConds),:);
%         LFPbySweep{2}=L;
%     end
% end
% 
% if size(LFPbySweep{1},1)~=length(unique(dLGNpsth.unitTrials{1}))
%     disp('Mismatched number of trials in LFPbySweep and dLGNpsth');
%     return
% end
% 
% disp('2');
% if length(LFPbySweep)>1
%     hippoLFP{1}=LFPbySweep{2};
% else
%     hippoLFP{1}=LFPbySweep{1};
% end
% thetaDiff=getHippoTheta(hippoLFP,ledConds,stimConds,usel_noLED,usealls,LFP_Fs);
% figure(); 
% hist(nanmean(thetaDiff,2),30);
% title('Histogram of theta diff');
% 
% cxSynch=[];
% if length(LFPbySweep)>1
%     cxLFP{1}=LFPbySweep{1};
%     cxSynch=getCortexSynch(cxLFP,ledConds,stimConds,usel_noLED,usealls,LFP_Fs);
% end
% 
% if length(LFPbySweep)>1
%     figure(); 
%     scatter(nanmean(thetaDiff,2),nanmean(cxSynch,2)); 
% end
% save([outputDir '\' 'thetaDiff_noLED.mat'],'thetaDiff');
% save([outputDir '\' 'cxSynch_noLED.mat'],'cxSynch');
% save([outputDir '\' 'dLGNpsth.mat'],'dLGNpsth');
% 
% thetaDiff=getHippoTheta(hippoLFP,ledConds,stimConds,[usel_noLED usel_LED],usealls,LFP_Fs);
% save([outputDir '\' 'thetaDiff.mat'],'thetaDiff');
% 
% thetaThresh=input('Theta difference threshold to use?');
% save([outputDir '\' 'thetaThresh.mat'],'thetaThresh');
% 
% noThetaTrials=nanmean(thetaDiff,2)<thetaThresh;
save([outputDir '\' 'noThetaTrials.mat'],'noThetaTrials');

disp('3');
tri=dLGNpsth.unitTrials{1};
if length(noThetaTrials)~=length(tri)
    disp('Length of noThetaTrials does not match number of trials in dLGNpsth');
    return
end
filt_psth=filtPSTH(dLGNpsth,tri(noThetaTrials'));
filt_psth_theta=filtPSTH(dLGNpsth,tri(noThetaTrials'~=1));

% No LED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usel=usel_noLED;

% uses=   {{1 1 1}; {1 1}};
% uses_tri={{[13:72]-12 [141:212]-12 [277:316]-12}; {[73:140]-12 [213:276]-12}};

for i=1:length(uses)
    currs=uses{i};
    currtri=uses_tri{i};
    allcurr_psth=[];
    allcurr_psth_theta=[];
    allcurr_psth_alltheta=[];
    for j=1:length(currs)
        subcurrs=currs{j};
        subcurrtri=currtri{j};
        
        sforsub=filt_psth.unitStimcond{1};
        triforsub=filt_psth.unitTrials{1};
        curr_psth=filtPSTH(filt_psth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        sforsub=filt_psth_theta.unitStimcond{1};
        triforsub=filt_psth_theta.unitTrials{1};
        curr_psth_theta=filtPSTH(filt_psth_theta,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        sforsub=dLGNpsth.unitStimcond{1};
        triforsub=dLGNpsth.unitTrials{1};
        curr_psth_alltheta=filtPSTH(dLGNpsth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
   
        allcurr_psth=concatPSTHs(allcurr_psth,curr_psth);
        allcurr_psth_theta=concatPSTHs(allcurr_psth_theta,curr_psth_theta);
        allcurr_psth_alltheta=concatPSTHs(allcurr_psth_alltheta,curr_psth_alltheta);
    end    
    [noTheta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),0,usealls,usel);
    [theta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel);
    [alltheta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth_alltheta,ones(length(dLGNpsth.psths),1),zeros(length(dLGNpsth.psths),1),0,usealls,usel);
end
save([outputDir '\' 'noTheta_trialAv_temp_noLED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_noLED.mat'],'theta_trialAv_temp');
save([outputDir '\' 'alltheta_trialAv_temp_noLED.mat'],'alltheta_trialAv_temp');
% noTheta_trialAv=avResponses(noTheta_trialAv_temp);
% theta_trialAv=avResponses(theta_trialAv_temp);
% save([outputDir '\' 'noTheta_trialAv_noLED.mat'],'noTheta_trialAv');
% save([outputDir '\' 'theta_trialAv_noLED.mat'],'theta_trialAv');

% [noTheta.allS]=sub_getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel);
% save([outputDir '\' 'noTheta_noLED.mat'],'noTheta');
% [theta.allS]=sub_getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel);
% save([outputDir '\' 'theta_noLED.mat'],'theta');

disp('4');
% With LED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usel=usel_LED;
for i=1:length(uses)
    currs=uses{i};
    currtri=uses_tri{i};
    allcurr_psth=[];
    allcurr_psth_theta=[];
    allcurr_psth_alltheta=[];
    for j=1:length(currs)
        subcurrs=currs{j};
        subcurrtri=currtri{j};
        
        sforsub=filt_psth.unitStimcond{1};
        triforsub=filt_psth.unitTrials{1};
        curr_psth=filtPSTH(filt_psth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        sforsub=filt_psth_theta.unitStimcond{1};
        triforsub=filt_psth_theta.unitTrials{1};
        curr_psth_theta=filtPSTH(filt_psth_theta,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        sforsub=dLGNpsth.unitStimcond{1};
        triforsub=dLGNpsth.unitTrials{1};
        curr_psth_alltheta=filtPSTH(dLGNpsth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));  
        
        allcurr_psth=concatPSTHs(allcurr_psth,curr_psth);
        allcurr_psth_theta=concatPSTHs(allcurr_psth_theta,curr_psth_theta);
        allcurr_psth_alltheta=concatPSTHs(allcurr_psth_alltheta,curr_psth_alltheta);
    end    
    [noTheta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),0,usealls,usel);
    [theta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth_theta,ones(length(allcurr_psth_theta.psths),1),zeros(length(allcurr_psth_theta.psths),1),0,usealls,usel);
    [alltheta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth_alltheta,ones(length(dLGNpsth.psths),1),zeros(length(dLGNpsth.psths),1),0,usealls,usel);
end
save([outputDir '\' 'noTheta_trialAv_temp_LED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_LED.mat'],'theta_trialAv_temp');
save([outputDir '\' 'alltheta_trialAv_temp_LED.mat'],'alltheta_trialAv_temp');
% noTheta_trialAv=avResponses(noTheta_trialAv_temp);
% theta_trialAv=avResponses(theta_trialAv_temp);
% save([outputDir '\' 'noTheta_trialAv_LED.mat'],'noTheta_trialAv');
% save([outputDir '\' 'theta_trialAv_LED.mat'],'theta_trialAv');

% [noTheta.allS]=sub_getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel);
% save([outputDir '\' 'noTheta_LED.mat'],'noTheta');
% [theta.allS]=sub_getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel);
% save([outputDir '\' 'theta_LED.mat'],'theta');

end

function outresponse=avResponses(response)

HFa=zeros(size(response(1).HFa));
LFa=zeros(size(response(1).LFa));
F1amp=zeros(size(response(1).F1amp));
allpower=zeros(size(response(1).allpower));
for i=1:length(response)
    if ~isnan(response(i).HFa)
        HFa=HFa+response(i).HFa;
    end
    if ~isnan(response(i).LFa)
        LFa=LFa+response(i).LFa;
    end
    if ~isnan(response(i).F1amp)
        F1amp=F1amp+response(i).F1amp;
    end
    if ~isnan(response(i).allpower)
        allpower=allpower+response(i).allpower;
    end
end
HFa=HFa./length(response);
LFa=LFa./length(response);
F1amp=F1amp./length(response);
allpower=allpower./length(response);
outresponse.HFa=HFa;
outresponse.LFa=LFa;
outresponse.F1amp=F1amp;
outresponse.allpower=allpower;

end


function newPSTH=concatPSTHs(struct1,struct2)

if isempty(struct1)
    newPSTH=struct2;
    return
elseif isempty(struct2)
    newPSTH=struct1;
    return
else
    if length(struct1.t)~=length(struct2.t)
        disp('t in psths does not match');
        return
    end
    newPSTH.t=struct1.t;
    newPSTH.psths=cell(length(struct1.psths),1);
    for i=1:length(struct1.psths)
        newPSTH.psths{i}=[struct1.psths{i}; struct2.psths{i}];
        newPSTH.unitTrials{i}=[struct1.unitTrials{i} struct2.unitTrials{i}];
        newPSTH.unitStimcond{i}=[struct1.unitStimcond{i} struct2.unitStimcond{i}];
        newPSTH.unitLED{i}=[struct1.unitLED{i} struct2.unitLED{i}];
    end
end

end


function [allS]=sub_getHFandLFalphaResponses(psth1,likes1,sigCells,avBeforeSpec,uses,usel)

allpower=[];
allS=[];
HFa=[];
LFa=[];
F1amp=[];

% usel=[0];
% uses=[1];
% HFrange=[11.5 12.5];
HFrange=[11.5 20];
LFrange=[4 6];
F1range=[2.5 3.5];
allrange=[10 99.5];
% avBeforeSpec=1;

movingwin=[1 0.05];
params.tapers=[5 6];
params.Fs=1/(psth1.t(2)-psth1.t(1));
% params.fpass=[1 30];
params.fpass=[1 100];
params.trialave=0;

useUnits=find(likes1==1 & sigCells<=0.05);
% useUnits=find(likes1==1 & sigCells<=0.5);
allS.t=[];
allS.f=[];
allS.S=cell(1,length(useUnits));
for i=1:length(useUnits)
    p=psth1.psths{useUnits(i)};
    l=psth1.unitLED{useUnits(i)};
    s=psth1.unitStimcond{useUnits(i)};
    if avBeforeSpec==1
        p=nanmean(p(ismember(l,usel) & ismember(s,uses),:),1)';
        [S,t,f]=mtspecgrampb(p,movingwin,params);
    else
        [S,t,f]=mtspecgrampb(p(ismember(l,usel) & ismember(s,uses),:)',movingwin,params);
    end
    allS.S{i}=S;
    allS.t=t;
    allS.f=f;
%     HFa(i,:)=nanmean(S(:,f>=HFrange(1) & f<=HFrange(2)),2)';
%     LFa(i,:)=nanmean(S(:,f>=LFrange(1) & f<=LFrange(2)),2)';
%     F1amp(i,:)=nanmean(S(:,f>=F1range(1) & f<=F1range(2)),2)';
%     allpower(i,:)=nanmean(S(:,f>=allrange(1) & f<=allrange(2)),2)';
end

end

