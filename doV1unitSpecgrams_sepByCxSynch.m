function doV1unitSpecgrams_sepByCxSynch(outputDir,dLGNpsth,usel_noLED,usel_LED,uses,uses_tri,noThetaTrials)


% Mawake2
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}};
% uses_tri={{[1:168]}; {[1:168]}; {[1:168]}; {[1:168]}; {[1:168]}; {[1:168]}; {[1:168]}; {[1:168]}; {[1:168]}};
% Mawake3
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}};
% uses_tri={{[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}};
% Mawake6
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}};
% uses_tri={{[1:216 409:528]}; {[1:216 409:528]}; {[1:216 409:528]}; {[1:216 409:528]}; {[1:216 409:528]}; {[1:216 409:528]}; {[1:216 409:528]}; {[1:216 409:528]}; {[1:216 409:528]}};
% Mawake7
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}};
% uses_tri={{[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}};
% Mawake8
% uses=   {{1}; {2}; {3}};
% uses_tri={{[1:972]}; {[1:972]}; {[1:972]}};
% Mawake18
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}; {10}; {11}; {12}; {13}; {14}; {15}; {16}};
% uses_tri={{[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}};
% Mawake24
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}; {10}; {11}; {12}; {13}; {14}; {15}; {16}};
% uses_tri={{[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}};
% Mawake25
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}};
% uses_tri={{[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}; {[1:132]}};
% Mawake26
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}; {10}; {11}; {12}; {13}; {14}; {15}; {16}};
% uses_tri={{[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}; {[1:420]}};
% Mawake27
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}};
% uses_tri={{[1:348]}; {[1:348]}; {[1:348]}; {[1:348]}; {[1:348]}; {[1:348]}; {[1:348]}; {[1:348]}; {[1:348]}};
% Mawake35
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}};
% uses_tri={{[1:96]}; {[1:96]}; {[1:96]}; {[1:96]}; {[1:96]}; {[1:96]}; {[1:96]}; {[1:96]}; {[1:96]}};
% Mawake44
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}};
% uses_tri={{[1:244 805:1048]}; {[1:244 805:1048]}; {[1:244 805:1048]}; {[1:244 805:1048]}; {[1:244 805:1048]}; {[1:244 805:1048]}; {[1:244 805:1048]}; {[1:244 805:1048]}; {[1:244 805:1048]}};
% Mawake96
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}; {10}; {11}; {12}; {13}; {14}; {15}; {16}; {17}; {18}; {19}; {20}; {21}; {22}; {23}; {24}; {25}; {26}; {27}; {28}; {29}; {30}; {31}; {32}};
% uses_tri={{[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}; {[1:396]}};


%uses=unique(dLGNpsth.unitStimcond{1});
usealls=1:10000;

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
    [noTheta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),0,usealls,usel,0);
    [theta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel,0);
    [alltheta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth_alltheta,ones(length(dLGNpsth.psths),1),zeros(length(dLGNpsth.psths),1),0,usealls,usel,0);
end
save([outputDir '\' 'noTheta_trialAv_temp_noLED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_noLED.mat'],'theta_trialAv_temp');
save([outputDir '\' 'alltheta_trialAv_temp_noLED.mat'],'alltheta_trialAv_temp');
% noTheta_trialAv=avResponses(noTheta_trialAv_temp);
% theta_trialAv=avResponses(theta_trialAv_temp);
% save([outputDir '\' 'noTheta_trialAv_noLED.mat'],'noTheta_trialAv');
% save([outputDir '\' 'theta_trialAv_noLED.mat'],'theta_trialAv');

[noTheta.allS]=sub_getHFandLFalphaResponses(dLGNpsth,ones(length(dLGNpsth.psths),1),zeros(length(dLGNpsth.psths),1),0,usealls,usel,1);
save([outputDir '\' 'allTheta_noLED.mat'],'noTheta');
[theta.allS]=sub_getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel,1);
save([outputDir '\' 'theta_noLED.mat'],'theta');
[noTheta.allS]=sub_getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel,1);
save([outputDir '\' 'noTheta_noLED.mat'],'noTheta');

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
    [noTheta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),0,usealls,usel,0);
    [theta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth_theta,ones(length(allcurr_psth_theta.psths),1),zeros(length(allcurr_psth_theta.psths),1),0,usealls,usel,0);
    [alltheta_trialAv_temp(i).allS]=sub_getHFandLFalphaResponses(allcurr_psth_alltheta,ones(length(dLGNpsth.psths),1),zeros(length(dLGNpsth.psths),1),0,usealls,usel,0);
end
save([outputDir '\' 'noTheta_trialAv_temp_LED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_LED.mat'],'theta_trialAv_temp');
save([outputDir '\' 'alltheta_trialAv_temp_LED.mat'],'alltheta_trialAv_temp');
% noTheta_trialAv=avResponses(noTheta_trialAv_temp);
% theta_trialAv=avResponses(theta_trialAv_temp);
% save([outputDir '\' 'noTheta_trialAv_LED.mat'],'noTheta_trialAv');
% save([outputDir '\' 'theta_trialAv_LED.mat'],'theta_trialAv');

[noTheta.allS]=sub_getHFandLFalphaResponses(dLGNpsth,ones(length(dLGNpsth.psths),1),zeros(length(dLGNpsth.psths),1),0,usealls,usel,1);
save([outputDir '\' 'allTheta_LED.mat'],'noTheta');
[theta.allS]=sub_getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel,1);
save([outputDir '\' 'theta_LED.mat'],'theta');
[noTheta.allS]=sub_getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel,1);
save([outputDir '\' 'noTheta_LED.mat'],'noTheta');

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


function [allS]=sub_getHFandLFalphaResponses(psth1,likes1,sigCells,avBeforeSpec,uses,usel,trialave)

deleteOnset=0;

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
params.fpass=[1 50];
params.trialave=trialave;

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
        if deleteOnset==1
            p=nanmean(p(ismember(l,usel) & ismember(s,uses),allS.t<4.02 | allS.t>4.2),1)';
        else
            p=nanmean(p(ismember(l,usel) & ismember(s,uses),:),1)';
        end
        [S,t,f]=mtspecgrampb(p,movingwin,params);
    else
        if deleteOnset==1
            p=p(ismember(l,usel) & ismember(s,uses),allS.t<4.02 | allS.t>4.2)';
        else
            p=p(ismember(l,usel) & ismember(s,uses),:)';
        end
%         [S,t,f]=mtspecgrampb(p(ismember(l,usel) & ismember(s,uses),:)',movingwin,params);
        [S,t,f]=mtspecgrampb(p,movingwin,params);
    end
    allS.S{i}=S;
    if deleteOnset==1
        allS.t=t(allS.t<4.02 | allS.t>4.2);
    else
        allS.t=t;
    end
    allS.f=f;
%     HFa(i,:)=nanmean(S(:,f>=HFrange(1) & f<=HFrange(2)),2)';
%     LFa(i,:)=nanmean(S(:,f>=LFrange(1) & f<=LFrange(2)),2)';
%     F1amp(i,:)=nanmean(S(:,f>=F1range(1) & f<=F1range(2)),2)';
%     allpower(i,:)=nanmean(S(:,f>=allrange(1) & f<=allrange(2)),2)';
end

end

