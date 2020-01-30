function combine_dLGN_vs_ntsr1_analysis(datadir,addToName,dLGNdir,V1dir,workspace,sortingWindow)

placeStimOnsetAt=4; % seconds from trial onset
saveFieldByField=true;
% maxFreq=330; % in Hz
maxFreq_override=50; % in Hz
maxTime_override=12; % in seconds from start
normByIntegralFreqRange=[10 45]; % in Hz, will divide through by integral of spectrum over this window
alphaRange=[10 14];
% sortingWindow=[3 3.6]; % time window used to sort trials into low and high
% sortingWindow=[4 5]; % time window used to sort trials into low and high

if isempty(dLGNdir)
    dLGNdir='\dLGN\specs\';
end
if isempty(V1dir)
    V1dir='\V1 Ntsr1\specs\';
end

if isempty(workspace)
    if iscell(datadir)
        
        all_units_noTheta_lowF1_Ntsr1=[];
        all_units_noTheta_highF1_Ntsr1=[];
        all_units_theta_lowF1_Ntsr1=[];
        all_units_theta_highF1_Ntsr1=[];
        
        for i=1:length(datadir)
            d=datadir{i};
            a=load([d '\' 'stimWindow.mat']);
            stimWindow=a.stimWindow;
            
            a=load([d dLGNdir 'noTheta_lowF1_dLGN' addToName]);
            noTheta_lowF1_dLGN=a.takeLow;
            a=load([d dLGNdir 'noTheta_highF1_dLGN' addToName]);
            noTheta_highF1_dLGN=a.takeHigh;
            a=load([d dLGNdir 'theta_lowF1_dLGN' addToName]);
            theta_lowF1_dLGN=a.takeLow;
            a=load([d dLGNdir 'theta_highF1_dLGN' addToName]);
            theta_highF1_dLGN=a.takeHigh;
            if saveFieldByField==true
                a=loadVar([d dLGNdir],'noTheta_noLED',saveFieldByField,'noTheta');
                a.noTheta.allS.S=loadVar([d dLGNdir],'noTheta_noLED_allS_S',saveFieldByField,'S');
                b=loadVar([d dLGNdir],'noTheta_trialAv_temp_noLED',saveFieldByField,'noTheta_trialAv_temp');
                a.noTheta.allS.t=b.noTheta_trialAv_temp(1).allS.t;
                a.noTheta.allS.f=b.noTheta_trialAv_temp(1).allS.f;
                noTheta=a.noTheta;
            else
                a=load([d dLGNdir 'noTheta_noLED']);
                noTheta=a.noTheta;
            end
            
            a=load([d V1dir 'noTheta_lowF1_Ntsr1' addToName]);
            noTheta_lowF1_Ntsr1=a.takeLow;
            a=load([d V1dir 'noTheta_highF1_Ntsr1' addToName]);
            noTheta_highF1_Ntsr1=a.takeHigh;
            a=load([d V1dir 'theta_lowF1_Ntsr1' addToName]);
            theta_lowF1_Ntsr1=a.takeLow;
            a=load([d V1dir 'theta_highF1_Ntsr1' addToName]);
            theta_highF1_Ntsr1=a.takeHigh;
            if saveFieldByField==true
                a=loadVar([d V1dir],'noTheta_noLED',saveFieldByField,'noTheta');
                a.noTheta.allS.S=loadVar([d V1dir],'noTheta_noLED_allS_S',saveFieldByField,'S');
                b=loadVar([d V1dir],'noTheta_trialAv_temp_noLED',saveFieldByField,'noTheta_trialAv_temp');
                a.noTheta.allS.t=b.noTheta_trialAv_temp(1).allS.t;
                a.noTheta.allS.f=b.noTheta_trialAv_temp(1).allS.f;
                noTheta_V1=a.noTheta;
            else
                a=load([d V1dir 'noTheta_noLED']);
                noTheta_V1=a.noTheta;
            end
            
            a=load([d V1dir 'noTheta_lowF1_Ntsr1_individUnits' addToName]);
            noTheta_lowF1_Ntsr1_individUnits=a.takeLow;
            a=load([d V1dir 'noTheta_highF1_Ntsr1_individUnits' addToName]);
            noTheta_highF1_Ntsr1_individUnits=a.takeHigh;
            a=load([d V1dir 'theta_lowF1_Ntsr1_individUnits' addToName]);
            theta_lowF1_Ntsr1_individUnits=a.takeLow;
            a=load([d V1dir 'theta_highF1_Ntsr1_individUnits' addToName]);
            theta_highF1_Ntsr1_individUnits=a.takeHigh;
            
            noTheta_lowF1_dLGN=alignToStim(noTheta_lowF1_dLGN,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
            noTheta_highF1_dLGN=alignToStim(noTheta_highF1_dLGN,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
            theta_lowF1_dLGN=alignToStim(theta_lowF1_dLGN,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
            theta_highF1_dLGN=alignToStim(theta_highF1_dLGN,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
            
            noTheta_lowF1_Ntsr1=alignToStim(noTheta_lowF1_Ntsr1,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
            noTheta_highF1_Ntsr1=alignToStim(noTheta_highF1_Ntsr1,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
            theta_lowF1_Ntsr1=alignToStim(theta_lowF1_Ntsr1,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
            theta_highF1_Ntsr1=alignToStim(theta_highF1_Ntsr1,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
            
            if i==1
                all_noTheta_lowF1_dLGN=zeros(size(noTheta_lowF1_dLGN));
                all_noTheta_highF1_dLGN=zeros(size(noTheta_lowF1_dLGN));
                all_theta_lowF1_dLGN=zeros(size(noTheta_lowF1_dLGN));
                all_theta_highF1_dLGN=zeros(size(noTheta_lowF1_dLGN));
                
                all_noTheta_lowF1_Ntsr1=zeros(size(noTheta_lowF1_Ntsr1));
                all_noTheta_highF1_Ntsr1=zeros(size(noTheta_lowF1_Ntsr1));
                all_theta_lowF1_Ntsr1=zeros(size(noTheta_lowF1_Ntsr1));
                all_theta_highF1_Ntsr1=zeros(size(noTheta_lowF1_Ntsr1));
            end
            
            
            all_noTheta_lowF1_dLGN=addAll(all_noTheta_lowF1_dLGN,length(noTheta.allS.S)*noTheta_lowF1_dLGN);
            all_noTheta_highF1_dLGN=addAll(all_noTheta_highF1_dLGN,length(noTheta.allS.S)*noTheta_highF1_dLGN);
            all_theta_lowF1_dLGN=addAll(all_theta_lowF1_dLGN,length(noTheta.allS.S)*theta_lowF1_dLGN);
            all_theta_highF1_dLGN=addAll(all_theta_highF1_dLGN,length(noTheta.allS.S)*theta_highF1_dLGN);
            
            all_noTheta_lowF1_Ntsr1=addAll(all_noTheta_lowF1_Ntsr1,length(noTheta.allS.S)*noTheta_lowF1_Ntsr1);
            all_noTheta_highF1_Ntsr1=addAll(all_noTheta_highF1_Ntsr1,length(noTheta.allS.S)*noTheta_highF1_Ntsr1);
            all_theta_lowF1_Ntsr1=addAll(all_theta_lowF1_Ntsr1,length(noTheta.allS.S)*theta_lowF1_Ntsr1);
            all_theta_highF1_Ntsr1=addAll(all_theta_highF1_Ntsr1,length(noTheta.allS.S)*theta_highF1_Ntsr1);
            
            all_units_noTheta_lowF1_Ntsr1=cat(3,all_units_noTheta_lowF1_Ntsr1,noTheta_lowF1_Ntsr1_individUnits);
            all_units_noTheta_highF1_Ntsr1=cat(3,all_units_noTheta_highF1_Ntsr1,noTheta_highF1_Ntsr1_individUnits);
            all_units_theta_lowF1_Ntsr1=cat(3,all_units_theta_lowF1_Ntsr1,theta_lowF1_Ntsr1_individUnits);
            all_units_theta_highF1_Ntsr1=cat(3,all_units_theta_highF1_Ntsr1,theta_highF1_Ntsr1_individUnits);
            
        end
    else
        error('expected datadir to be a cell array');
    end
else
    load(workspace);
end

maxFreq=maxFreq_override;
maxTime=maxTime_override;
all_units_noTheta_lowF1_Ntsr1=backup_all_units_noTheta_lowF1_Ntsr1;
all_units_noTheta_highF1_Ntsr1=backup_all_units_noTheta_highF1_Ntsr1;
all_units_theta_lowF1_Ntsr1=backup_all_units_theta_lowF1_Ntsr1;
all_units_theta_highF1_Ntsr1=backup_all_units_theta_highF1_Ntsr1;

% sigMask_noTheta=calculateSignificanceMask(all_units_noTheta_lowF1_Ntsr1,all_units_noTheta_highF1_Ntsr1);
% sigMask_theta=calculateSignificanceMask(all_units_theta_lowF1_Ntsr1,all_units_theta_highF1_Ntsr1);

all_noTheta_lowF1_dLGN=all_noTheta_lowF1_dLGN(noTheta.allS.t<maxTime,noTheta.allS.f<=maxFreq);
all_noTheta_highF1_dLGN=all_noTheta_highF1_dLGN(noTheta.allS.t<maxTime,noTheta.allS.f<=maxFreq);
all_theta_lowF1_dLGN=all_theta_lowF1_dLGN(noTheta.allS.t<maxTime,noTheta.allS.f<=maxFreq);
all_theta_highF1_dLGN=all_theta_highF1_dLGN(noTheta.allS.t<maxTime,noTheta.allS.f<=maxFreq);

all_noTheta_lowF1_Ntsr1=all_noTheta_lowF1_Ntsr1(noTheta.allS.t<maxTime,noTheta_V1.allS.f<=maxFreq);
all_noTheta_highF1_Ntsr1=all_noTheta_highF1_Ntsr1(noTheta.allS.t<maxTime,noTheta_V1.allS.f<=maxFreq);
all_theta_lowF1_Ntsr1=all_theta_lowF1_Ntsr1(noTheta.allS.t<maxTime,noTheta_V1.allS.f<=maxFreq);
all_theta_highF1_Ntsr1=all_theta_highF1_Ntsr1(noTheta.allS.t<maxTime,noTheta_V1.allS.f<=maxFreq);

% sig_noTheta_lowF1_Ntsr1=addSigMask(all_noTheta_lowF1_Ntsr1,sigMask_noTheta,noTheta,0.8);
% sig_noTheta_highF1_Ntsr1=addSigMask(all_noTheta_highF1_Ntsr1,sigMask_noTheta,noTheta,0.8);
% sig_theta_lowF1_Ntsr1=addSigMask(all_theta_lowF1_Ntsr1,sigMask_theta,noTheta,0.8);
% sig_theta_highF1_Ntsr1=addSigMask(all_theta_highF1_Ntsr1,sigMask_theta,noTheta,0.8);

backup_all_units_noTheta_lowF1_Ntsr1=all_units_noTheta_lowF1_Ntsr1;
backup_all_units_noTheta_highF1_Ntsr1=all_units_noTheta_highF1_Ntsr1;
backup_all_units_theta_lowF1_Ntsr1=all_units_theta_lowF1_Ntsr1;
backup_all_units_theta_highF1_Ntsr1=all_units_theta_highF1_Ntsr1;
% Do whitening
% all_units_noTheta_lowF1_Ntsr1=whitenAll(all_units_noTheta_lowF1_Ntsr1,1,0);
% all_units_noTheta_highF1_Ntsr1=whitenAll(all_units_noTheta_highF1_Ntsr1,1,1);
% all_units_theta_lowF1_Ntsr1=whitenAll(all_units_theta_lowF1_Ntsr1,1,0);
% all_units_theta_highF1_Ntsr1=whitenAll(all_units_theta_highF1_Ntsr1,1,1);
% all_units_noTheta_lowF1_Ntsr1=whitenAll(all_units_noTheta_lowF1_Ntsr1,0,0);
% all_units_noTheta_highF1_Ntsr1=whitenAll(all_units_noTheta_highF1_Ntsr1,0,1);
% all_units_theta_lowF1_Ntsr1=whitenAll(all_units_theta_lowF1_Ntsr1,0,0);
% all_units_theta_highF1_Ntsr1=whitenAll(all_units_theta_highF1_Ntsr1,0,1);

% Norm
% all_units_noTheta_lowF1_Ntsr1=normAll(all_units_noTheta_lowF1_Ntsr1);
% all_units_noTheta_highF1_Ntsr1=normAll(all_units_noTheta_highF1_Ntsr1);
% all_units_theta_lowF1_Ntsr1=normAll(all_units_theta_lowF1_Ntsr1);
% all_units_theta_highF1_Ntsr1=normAll(all_units_theta_highF1_Ntsr1);

% Find units that have alpha power that correlates with Ntsr1
% i.e., alpha power is lower for "low F1" and higher for "high F1", if did
% sorting based on alpha power in Ntsr1 unit
% diffThresh=0.15; 
% diffThresh=0.13; 
% diffThresh=0.11; % good! when paired with diffThresh for vis evoked 0,n=9
% diffThresh=0.05; diffThresh=2;
diffThresh=0.05; 
versusRange=[20 40];
[useUnits_noTheta_alphaMatch,diffHighMinusLow_noTheta]=isAlphaHigherForHigh(all_units_noTheta_lowF1_Ntsr1,all_units_noTheta_highF1_Ntsr1,alphaRange,sortingWindow,noTheta_V1.allS.t,noTheta_V1.allS.f,diffThresh,versusRange);
[useUnits_theta_alphaMatch,diffHighMinusLow_theta]=isAlphaHigherForHigh(all_units_theta_lowF1_Ntsr1,all_units_theta_highF1_Ntsr1,alphaRange,sortingWindow,noTheta_V1.allS.t,noTheta_V1.allS.f,diffThresh,versusRange);
% useUnits_noTheta_alphaMatch=diffHighMinusLow_noTheta<0.3 & diffHighMinusLow_noTheta>0.15;
% useUnits_noTheta_alphaMatch=-diffHighMinusLow_noTheta<0.11;

diffThresh=0; 
[useUnits_vis_noTheta,diffEvMinusSpont_noTheta]=isVisuallyResponsive(all_units_theta_lowF1_Ntsr1,all_units_theta_highF1_Ntsr1,[3 3.6],[4 6.5],t,f,diffThresh,[1.5 2.5]);
[useUnits_vis_theta,diffEvMinusSpont_theta]=isVisuallyResponsive(all_units_theta_lowF1_Ntsr1,all_units_theta_highF1_Ntsr1,[3 3.6],[4 6.5],t,f,diffThresh,[1.5 2.5]);
% useUnits_vis_noTheta=diffEvMinusSpont_noTheta<diffThresh;
useUnits_vis_noTheta=diffEvMinusSpont_noTheta<=1.25 & diffEvMinusSpont_noTheta>0;
% useUnits_vis_noTheta=diffEvMinusSpont_noTheta>6;

[alphaDiff,F1responseForGLM,meanRateForGLM,visEvForGLM]=plotXAsFuncOfY(all_units_noTheta_lowF1_Ntsr1,all_units_noTheta_highF1_Ntsr1,alphaRange,versusRange,sortingWindow,noTheta,maxTime,[-0.8+0.01:0.02:0.8-0.01],[0 50],diffEvMinusSpont_noTheta);

[meanRateDiff,F1response]=plotXAsFuncOfY(all_units_noTheta_lowF1_Ntsr1,all_units_noTheta_highF1_Ntsr1,[0 50],[],sortingWindow,noTheta,maxTime,[-0.8+0.01:0.02:0.8-0.01],[],[]);

mdl=fitglm([alphaDiff meanRateForGLM visEvForGLM],F1responseForGLM)

[~,diffTest]=isAlphaHigherForHigh(all_units_noTheta_lowF1_Ntsr1,all_units_noTheta_highF1_Ntsr1,alphaRange,sortingWindow,noTheta_V1.allS.t,noTheta_V1.allS.f,diffThresh,versusRange);
[~,diffTest2]=isAlphaHigherForHigh(all_units_noTheta_lowF1_Ntsr1,all_units_noTheta_highF1_Ntsr1,[2.5 3.5],sortingWindow,noTheta_V1.allS.t,noTheta_V1.allS.f,diffThresh,[]);
figure(); scatter(diffTest2,diffTest); xlabel('low frequency mean rate'); ylabel('relative alpha');
[r,p]=corrcoef(diffTest2(~isnan(diffTest2) & ~isnan(diffTest)),diffTest(~isnan(diffTest2) & ~isnan(diffTest)))



useUnits_noTheta=useUnits_noTheta_alphaMatch & useUnits_vis_noTheta;
useUnits_theta=useUnits_theta_alphaMatch & useUnits_vis_theta;

disp('n units no theta');
disp(nansum(useUnits_noTheta));

t=noTheta_V1.allS.t;
f=noTheta_V1.allS.f;
% Plot high vs low Ntsr1 for no theta, use units only
plotSpecgramsComparingHighAndLow(all_units_noTheta_lowF1_Ntsr1,all_units_noTheta_highF1_Ntsr1,useUnits_noTheta,t,f,maxTime,maxFreq,'no theta Ntsr1',false);
% Plot high vs low Ntsr1 for theta, use units only
plotSpecgramsComparingHighAndLow(all_units_theta_lowF1_Ntsr1,all_units_theta_highF1_Ntsr1,useUnits_theta,t,f,maxTime,maxFreq,'theta Ntsr1',false);

powerInBand_low=takeFrequencyBand(all_units_noTheta_lowF1_Ntsr1(t<=maxTime,:,useUnits_noTheta),[2.5 3.5],noTheta_V1);
powerInBand_high=takeFrequencyBand(all_units_noTheta_highF1_Ntsr1(t<=maxTime,:,useUnits_noTheta),[2.5 3.5],noTheta_V1);
plotWStderr(t(t<=maxTime),powerInBand_low,powerInBand_high,'k','r',0); 
% plotWStderr(t(t<=maxTime),powerInBand_low./repmat(nansum(powerInBand_low,2),1,size(powerInBand_low,2)),powerInBand_high./repmat(nansum(powerInBand_high,2),1,size(powerInBand_high,2)),'k','r',0); % normalize each unit

% spectrum_noTheta_lowF1_Ntsr1=takeTimeSlice(all_units_noTheta_lowF1_Ntsr1,[3 3.7],noTheta_V1);
% spectrum_noTheta_highF1_Ntsr1=takeTimeSlice(all_units_noTheta_highF1_Ntsr1,[3 3.7],noTheta_V1);
% spectrum_theta_lowF1_Ntsr1=takeTimeSlice(all_units_theta_lowF1_Ntsr1,[3 3.7],noTheta_V1);
% spectrum_theta_highF1_Ntsr1=takeTimeSlice(all_units_theta_highF1_Ntsr1,[3 3.7],noTheta_V1);
% temp1=spectrum_noTheta_lowF1_Ntsr1./repmat(nansum(spectrum_noTheta_lowF1_Ntsr1(:,f>=normByIntegralFreqRange(1) & f<=normByIntegralFreqRange(2)),2),1,size(spectrum_noTheta_lowF1_Ntsr1,2));
% temp2=spectrum_noTheta_highF1_Ntsr1./repmat(nansum(spectrum_noTheta_highF1_Ntsr1(:,f>=normByIntegralFreqRange(1) & f<=normByIntegralFreqRange(2)),2),1,size(spectrum_noTheta_highF1_Ntsr1,2));
% signrank(nanmean(temp1(:,f>17 & f<22),2),nanmean(temp2(:,f>17 & f<22),2))
% nanmean(nanmean(temp1(:,f>17 & f<22),2))
% nanmean(nanmean(temp2(:,f>17 & f<22),2))
% y1=(spectrum_noTheta_lowF1_Ntsr1.*repmat((0.008/45)*f+0.026,size(spectrum_noTheta_lowF1_Ntsr1,1),1))./repmat(nansum(spectrum_noTheta_lowF1_Ntsr1(:,f>=normByIntegralFreqRange(1) & f<=normByIntegralFreqRange(2)),2),1,size(spectrum_noTheta_lowF1_Ntsr1,2));
% y2=(spectrum_noTheta_highF1_Ntsr1.*repmat((0.008/45)*f+0.026,size(spectrum_noTheta_lowF1_Ntsr1,1),1))./repmat(nansum(spectrum_noTheta_highF1_Ntsr1(:,f>=normByIntegralFreqRange(1) & f<=normByIntegralFreqRange(2)),2),1,size(spectrum_noTheta_lowF1_Ntsr1,2));
% y1(:,4)=y1(:,4)+0.00026;
% y1(:,1:3)=y1(:,1:3)+0.00037;
% y1(:,3)=y1(:,3)+0.00002;
% y2(:,4)=y2(:,4)+0.00026;
% y2(:,1:3)=y2(:,1:3)+0.00037;
% y2(:,3)=y2(:,3)+0.00002;
% plotWStderr(f,y1,y2,'k','r',0);
% title('no theta low vs high F1 de-trend');


temp1=whitenAll(all_units_noTheta_lowF1_Ntsr1,0,0);
temp2=whitenAll(all_units_noTheta_highF1_Ntsr1,0,1);
plotSpecgramsComparingHighAndLow(temp1,temp2,useUnits_noTheta,t,f,maxTime,maxFreq,'no theta Ntsr1',false);

temp1=whitenAll(all_units_theta_lowF1_Ntsr1,0,0);
temp2=whitenAll(all_units_theta_highF1_Ntsr1,0,1);
plotSpecgramsComparingHighAndLow(temp1,temp2,useUnits_theta,t,f,maxTime,maxFreq,'theta Ntsr1',false);

return

plotWStderr(f,spectrum_noTheta_lowF1_Ntsr1,spectrum_noTheta_highF1_Ntsr1,'k','r',0);
title('Before stim frequency spectrum noTheta low vs high F1');
plotWStderr(f,spectrum_theta_lowF1_Ntsr1,spectrum_theta_highF1_Ntsr1,'k','r',0);
title('Before stim frequency spectrum theta low vs high F1');

alphaBand_noTheta_lowF1_Ntsr1=takeFrequencyBand(all_units_noTheta_lowF1_Ntsr1,[17 22],noTheta_V1);
alphaBand_noTheta_highF1_Ntsr1=takeFrequencyBand(all_units_noTheta_highF1_Ntsr1,[17 22],noTheta_V1);
temp1=alphaBand_noTheta_lowF1_Ntsr1;
temp2=alphaBand_noTheta_highF1_Ntsr1;
% plotWStderr(t,temp1-repmat(nanmean(temp1(:,t<1),2),1,271),temp2-repmat(nanmean(temp2(:,t<1),2),1,271),'k','r',0);
plotWStderr(t,(temp1-repmat(nanmean(temp1(:,t<1),2),1,271))/nanmean(nanmean(spectrum_noTheta_lowF1_Ntsr1(:,f>=10 & f<=45),2)),(temp2-repmat(nanmean(temp2(:,t<1),2),1,271))/nanmean(nanmean(spectrum_noTheta_highF1_Ntsr1(:,f>=10 & f<=45),2)),'k','r',0);

alphaBand_noTheta_lowF1_Ntsr1=takeFrequencyBand(all_units_noTheta_lowF1_Ntsr1,[7 22],noTheta_V1);
alphaBand_noTheta_highF1_Ntsr1=takeFrequencyBand(all_units_noTheta_highF1_Ntsr1,[7 22],noTheta_V1);
alphaBand_theta_lowF1_Ntsr1=takeFrequencyBand(all_units_theta_lowF1_Ntsr1,[7 22],noTheta_V1);
alphaBand_theta_highF1_Ntsr1=takeFrequencyBand(all_units_theta_highF1_Ntsr1,[7 22],noTheta_V1);

plotWStderr(noTheta_V1.allS.t,alphaBand_noTheta_lowF1_Ntsr1,alphaBand_noTheta_highF1_Ntsr1,'k','r',0);
title('Alpha band noTheta low vs high F1');
plotWStderr(noTheta_V1.allS.t,alphaBand_theta_lowF1_Ntsr1,alphaBand_theta_highF1_Ntsr1,'k','r',0);
title('Alpha band theta low vs high F1');



% alphaBand_noTheta_lowF1_Ntsr1=takeFrequencyBand(all_units_noTheta_lowF1_Ntsr1,[7 22],noTheta);
% alphaBand_noTheta_highF1_Ntsr1=takeFrequencyBand(all_units_noTheta_highF1_Ntsr1,[7 22],noTheta);
% tempLow=alphaBand_noTheta_lowF1_Ntsr1; tempHigh=alphaBand_noTheta_highF1_Ntsr1;
% tempLow=tempLow-repmat(nanmin(tempLow(:,noTheta.allS.t<=2.7),[],2),1,271); tempLow=tempLow./repmat(nanmean(tempLow(:,noTheta.allS.t<12),2),1,271); 
% tempHigh=tempHigh-repmat(nanmin(tempHigh(:,noTheta.allS.t<=2.7),[],2),1,271); tempHigh=tempHigh./repmat(nanmean(tempHigh(:,noTheta.allS.t<12),2),1,271); 
% plotWStderr(noTheta.allS.t,tempLow,tempHigh,'k','r',0);
% signrank(nanmean(tempLow(:,noTheta.allS.t>=0 & noTheta.allS.t<=3.5),2),nanmean(tempHigh(:,noTheta.allS.t>=0 & noTheta.allS.t<=3.5),2))

backupall_noTheta_lowF1_dLGN=all_noTheta_lowF1_dLGN;
backupall_noTheta_highF1_dLGN=all_noTheta_highF1_dLGN;
backupall_theta_lowF1_dLGN=all_theta_lowF1_dLGN;
backupall_theta_highF1_dLGN=all_theta_highF1_dLGN;

backupall_noTheta_lowF1_Ntsr1=all_noTheta_lowF1_Ntsr1;
backupall_noTheta_highF1_Ntsr1=all_noTheta_highF1_Ntsr1;
backupall_theta_lowF1_Ntsr1=all_theta_lowF1_Ntsr1;
backupall_theta_highF1_Ntsr1=all_theta_highF1_Ntsr1;

all_noTheta_lowF1_dLGN=whitenAndNormSpecgram(all_noTheta_lowF1_dLGN,false,0);
all_noTheta_highF1_dLGN=whitenAndNormSpecgram(all_noTheta_highF1_dLGN,false,1);
all_theta_lowF1_dLGN=whitenAndNormSpecgram(all_theta_lowF1_dLGN,false,0);
all_theta_highF1_dLGN=whitenAndNormSpecgram(all_theta_highF1_dLGN,false,1);

all_noTheta_lowF1_Ntsr1=whitenAndNormSpecgram(all_noTheta_lowF1_Ntsr1,true,0);
all_noTheta_highF1_Ntsr1=whitenAndNormSpecgram(all_noTheta_highF1_Ntsr1,true,1);
all_theta_lowF1_Ntsr1=whitenAndNormSpecgram(all_theta_lowF1_Ntsr1,true,0);
all_theta_highF1_Ntsr1=whitenAndNormSpecgram(all_theta_highF1_Ntsr1,true,1);


figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_noTheta_lowF1_dLGN'); title('no theta dLGN low F1');
figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_noTheta_highF1_dLGN'); title('no theta dLGN high F1');
figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_theta_lowF1_dLGN'); title('theta dLGN low F1');
figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_theta_highF1_dLGN'); title('theta dLGN high F1');

figure(); imagesc(noTheta_V1.allS.t(noTheta_V1.allS.t<12),noTheta_V1.allS.f(noTheta_V1.allS.f<=50),all_noTheta_lowF1_Ntsr1'); title('no theta Ntsr1 low F1');
figure(); imagesc(noTheta_V1.allS.t(noTheta_V1.allS.t<12),noTheta_V1.allS.f(noTheta_V1.allS.f<=50),all_noTheta_highF1_Ntsr1'); title('no theta Ntsr1 high F1');
figure(); imagesc(noTheta_V1.allS.t(noTheta_V1.allS.t<12),noTheta_V1.allS.f(noTheta_V1.allS.f<=50),all_theta_lowF1_Ntsr1'); title('theta Ntsr1 low F1');
figure(); imagesc(noTheta_V1.allS.t(noTheta_V1.allS.t<12),noTheta_V1.allS.f(noTheta_V1.allS.f<=50),all_theta_highF1_Ntsr1'); title('theta Ntsr1 high F1');

% figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),sig_noTheta_lowF1_Ntsr1'); title('no theta Ntsr1 low F1 w sig mask');
% figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),sig_noTheta_highF1_Ntsr1'); title('no theta Ntsr1 high F1');
% figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),sig_theta_lowF1_Ntsr1'); title('theta Ntsr1 low F1');
% figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),sig_theta_highF1_Ntsr1'); title('theta Ntsr1 high F1');

end

function [concatX,concatY,concatThird,concatParam]=plotXAsFuncOfY(all_units_low,all_units_high,useFreqRange,versusRange,useTimeRange,noTheta,maxTime,tryThresh,returnThirdVarRange,param)

useUnits_noTheta_alphaMatch=cell(1,length(tryThresh)-1);
for i=1:length(tryThresh)-1
    diffThresh=tryThresh(i);
    [~,diffHighMinusLow_noTheta]=isAlphaHigherForHigh(all_units_low,all_units_high,useFreqRange,useTimeRange,noTheta.allS.t,noTheta.allS.f,diffThresh,versusRange);
    useUnits_noTheta_alphaMatch{i}=diffHighMinusLow_noTheta>=tryThresh(i) & diffHighMinusLow_noTheta<=tryThresh(i+1);
end
   
t=noTheta.allS.t;
powLow=cell(1,length(tryThresh)-1);
powHigh=cell(1,length(tryThresh)-1);
powThird=cell(1,length(tryThresh)-1);
highMinusLow=cell(1,length(tryThresh)-1);
powParam=cell(1,length(tryThresh)-1);
for i=1:length(tryThresh)-1
    useUnits_noTheta=useUnits_noTheta_alphaMatch{i};
    powerInBand_low=takeFrequencyBand(all_units_low(t<=maxTime,:,useUnits_noTheta),[2.5 3.5],noTheta);
    powLow{i}=nanmean(powerInBand_low(:,t>=4 & t<=6.5),2);
    powerInBand_high=takeFrequencyBand(all_units_high(t<=maxTime,:,useUnits_noTheta),[2.5 3.5],noTheta);
    powHigh{i}=nanmean(powerInBand_high(:,t>=4 & t<=6.5),2);
    highMinusLow{i}=powHigh{i}-powLow{i};
    if ~isempty(returnThirdVarRange)
        powerInBand_third_low=takeFrequencyBand(all_units_low(t<=maxTime,:,useUnits_noTheta),returnThirdVarRange,noTheta);
        powerInBand_third_high=takeFrequencyBand(all_units_high(t<=maxTime,:,useUnits_noTheta),returnThirdVarRange,noTheta);
        powThird{i}=nanmean(powerInBand_third_high(:,t>=4 & t<=6.5),2)-nanmean(powerInBand_third_low(:,t>=4 & t<=6.5),2);
    end
    if ~isempty(param)
        powParam{i}=param(useUnits_noTheta)';
    end
end

figure();
concatX=[];
concatY=[];
concatThird=[];
concatParam=[];
isInParamRange=[];
for i=1:length(tryThresh)-1
    if ~isempty(param)
        paramThresh=0;
        paramThresh2=1.25;
        paramTemp=powParam{i};
        takeParaming=paramTemp>paramThresh & paramTemp<=paramThresh2;
        x=ones(size(highMinusLow{i})).*nanmean([tryThresh(i) tryThresh(i+1)],2);
        y=highMinusLow{i};
        scatter(x(~takeParaming),y(~takeParaming),[],'k'); hold all;
        scatter(x(takeParaming),y(takeParaming),[],'r');
        scatter(x(paramTemp>paramThresh2),y(paramTemp>paramThresh2),[],'g');
    else
        scatter(ones(size(highMinusLow{i})).*nanmean([tryThresh(i) tryThresh(i+1)],2),highMinusLow{i});
    end
    concatX=[concatX; ones(size(highMinusLow{i})).*nanmean([tryThresh(i) tryThresh(i+1)],2)];
    concatY=[concatY; highMinusLow{i}];
    concatThird=[concatThird; powThird{i}];
    if ~isempty(param)
        concatParam=[concatParam; powParam{i}];
        isInParamRange=[isInParamRange; takeParaming];
    end
    hold on;
end

[r,p]=corrcoef(concatX,concatY)

[n,x]=hist(concatY(isInParamRange==1),100);
figure(); plot(x,n./nansum(n),'Color','r');
hold on;
[n,x]=hist(concatY(isInParamRange==0),300);
plot(x,n./nansum(n),'Color','k');

end

function [useUnits,diffEvMinusSpont]=isVisuallyResponsive(all_units_low,all_units_high,spontWindow,evWindow,t,f,diffThresh,visRange)

if length(size(all_units_low))>2
    spont=0.5*reshape(nanmean(nanmean(all_units_low(t>=spontWindow(1) & t<=spontWindow(2),f>=visRange(1) & f<=visRange(2),:),2),1),1,size(all_units_low,3)) + 0.5*reshape(nanmean(nanmean(all_units_high(t>=spontWindow(1) & t<=spontWindow(2),f>=visRange(1) & f<=visRange(2),:),2),1),1,size(all_units_high,3));
    ev=0.5*reshape(nanmean(nanmean(all_units_low(t>=evWindow(1) & t<=evWindow(2),f>=visRange(1) & f<=visRange(2),:),2),1),1,size(all_units_low,3)) + 0.5*reshape(nanmean(nanmean(all_units_high(t>=evWindow(1) & t<=evWindow(2),f>=visRange(1) & f<=visRange(2),:),2),1),1,size(all_units_high,3));
else
    spont=0.5*reshape(nanmean(nanmean(all_units_low(t>=spontWindow(1) & t<=spontWindow(2),f>=visRange(1) & f<=visRange(2)),2),1),1,1) + 0.5*reshape(nanmean(nanmean(all_units_high(t>=spontWindow(1) & t<=spontWindow(2),f>=visRange(1) & f<=visRange(2)),2),1),1,1);
    ev=0.5*reshape(nanmean(nanmean(all_units_low(t>=evWindow(1) & t<=evWindow(2),f>=visRange(1) & f<=visRange(2)),2),1),1,1) + 0.5*reshape(nanmean(nanmean(all_units_high(t>=evWindow(1) & t<=evWindow(2),f>=visRange(1) & f<=visRange(2)),2),1),1,1);
end
if isempty(diffThresh)
    useUnits=spont<ev;
    diffEvMinusSpont=ev-spont;
else
    diffEvMinusSpont=ev-spont;
    useUnits=diffEvMinusSpont>diffThresh;
end


end

function plotSpecgramsComparingHighAndLow(all_units_noTheta_lowF1_Ntsr1,all_units_noTheta_highF1_Ntsr1,useUnits_noTheta,t,f,maxTime,maxFreq,titAdd,peakNorm)

smoothAndInterp=true;

test=reshape(nanmean(all_units_noTheta_lowF1_Ntsr1(:,:,useUnits_noTheta),3),size(all_units_noTheta_lowF1_Ntsr1,1),size(all_units_noTheta_lowF1_Ntsr1,2));
test_high=reshape(nanmean(all_units_noTheta_highF1_Ntsr1(:,:,useUnits_noTheta),3),size(all_units_noTheta_highF1_Ntsr1,1),size(all_units_noTheta_highF1_Ntsr1,2));

if peakNorm==true
    test=peakNormSpecgram(test);
    test_high=peakNormSpecgram(test_high);
end

if smoothAndInterp==true
    temp=test_high(t<maxTime,f<=maxFreq)-test(t<maxTime,f<=maxFreq);
    temp=peakNormSpecgram(temp);
    smoothMat=smoothAndInterpSpecgram(temp');
    figure(); imagesc(t(t<maxTime),f(f<=maxFreq),smoothMat); title(['high minus low ' titAdd]);
    smoothMat1=smoothAndInterpSpecgram(test_high(t<maxTime,f<=maxFreq));
    smoothMat2=smoothAndInterpSpecgram(test(t<maxTime,f<=maxFreq));
    smoothMat=[smoothMat1 smoothMat2]';
    figure(); imagesc(t(t<maxTime),f(f<=maxFreq),smoothMat); title(['high then low ' titAdd]);
    smoothMat=smoothAndInterpSpecgram([test_high(t<maxTime,f<=maxFreq)]');
    figure(); imagesc(t(t<maxTime),f(f<=maxFreq),smoothMat); title(['high ' titAdd]);
    smoothMat=smoothAndInterpSpecgram([test(t<maxTime,f<=maxFreq)]');
    figure(); imagesc(t(t<maxTime),f(f<=maxFreq),smoothMat); title(['low ' titAdd]);
else
    figure(); imagesc(t(t<maxTime),f(f<=maxFreq),test_high(t<maxTime,f<=maxFreq)'-test(t<maxTime,f<=maxFreq)'); title(['high minus low ' titAdd]);
    figure(); imagesc(t(t<maxTime),f(f<=maxFreq),[test_high(t<maxTime,f<=maxFreq) test(t<maxTime,f<=maxFreq)]'); title(['high then low ' titAdd]);
    figure(); imagesc(t(t<maxTime),f(f<=maxFreq),[test_high(t<maxTime,f<=maxFreq)]'); title(['high ' titAdd]);
    figure(); imagesc(t(t<maxTime),f(f<=maxFreq),[test(t<maxTime,f<=maxFreq)]'); title(['low ' titAdd]);
end

end

function smoothMat=smoothAndInterpSpecgram(specgram)

smoothBy=3; 
K=ones(smoothBy);
temp=[specgram(1,:); specgram; specgram(end,:)];
temp=[temp(:,1) temp temp(:,end)];
smoothMat=conv2(temp,K,'same'); 
% smoothMat=smoothMat(smoothBy+1:end-smoothBy-1,smoothBy+1:end-smoothBy-1);
smoothMat=interp2(smoothMat,4);
% smoothBy=12;
smoothBy=22;
smoothMat=smoothMat(smoothBy+1:end-smoothBy-1,smoothBy+1:end-smoothBy-1);

end

function temp=peakNormSpecgram(temp)

temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));

end

function [useUnits,diffHighMinusLow]=isAlphaHigherForHigh(all_units_low,all_units_high,alphaRange,timeRange,t,f,diffThresh,versusRange)

if length(size(all_units_low))>2
    alphaForLow=reshape(nanmean(nanmean(all_units_low(t>=timeRange(1) & t<=timeRange(2),f>=alphaRange(1) & f<=alphaRange(2),:),2),1),1,size(all_units_low,3));
    alphaForHigh=reshape(nanmean(nanmean(all_units_high(t>=timeRange(1) & t<=timeRange(2),f>=alphaRange(1) & f<=alphaRange(2),:),2),1),1,size(all_units_high,3));
else
    alphaForLow=reshape(nanmean(nanmean(all_units_low(t>=timeRange(1) & t<=timeRange(2),f>=alphaRange(1) & f<=alphaRange(2)),2),1),1,1);
    alphaForHigh=reshape(nanmean(nanmean(all_units_high(t>=timeRange(1) & t<=timeRange(2),f>=alphaRange(1) & f<=alphaRange(2)),2),1),1,1);
end
if ~isempty(versusRange)
    if length(size(all_units_low))>2
        versusForLow=reshape(nanmean(nanmean(all_units_low(t>=timeRange(1) & t<=timeRange(2),f>=versusRange(1) & f<=versusRange(2),:),2),1),1,size(all_units_low,3));
        versusForHigh=reshape(nanmean(nanmean(all_units_high(t>=timeRange(1) & t<=timeRange(2),f>=versusRange(1) & f<=versusRange(2),:),2),1),1,size(all_units_high,3));
    else
        versusForLow=reshape(nanmean(nanmean(all_units_low(t>=timeRange(1) & t<=timeRange(2),f>=versusRange(1) & f<=versusRange(2)),2),1),1,1);
        versusForHigh=reshape(nanmean(nanmean(all_units_high(t>=timeRange(1) & t<=timeRange(2),f>=versusRange(1) & f<=versusRange(2)),2),1),1,1);
    end
    alphaForLow=alphaForLow./versusForLow;
    alphaForHigh=alphaForHigh./versusForHigh;
end
if isempty(diffThresh)
    useUnits=alphaForLow<alphaForHigh;
    diffHighMinusLow=alphaForHigh-alphaForLow;
else
    diffHighMinusLow=alphaForHigh-alphaForLow;
    useUnits=diffHighMinusLow>diffThresh;
end
    
end

function [y1,y2]=plotWStderr(x,y1,y2,c1,c2,doFill)

baseSub=false;
baseWindow=[1.485 3.535];

if baseSub==true
    y1=y1-repmat(nanmean(y1(:,x>=baseWindow(1) & x<=baseWindow(2)),2),1,size(y1,2));
    if ~isempty(y2)
        y2=y2-repmat(nanmean(y2(:,x>=baseWindow(1) & x<=baseWindow(2)),2),1,size(y2,2));
    end
end

figure();
plot(x,nanmean(y1,1),'Color',c1);
hold on;
if doFill==1
    fill([x fliplr(x)],[nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)) fliplr(nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)))],[0.5 0.5 0.5]);
end
plot(x,nanmean(y1,1),'Color',c1);
plot(x,nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);
plot(x,nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);

if ~isempty(y2)
    plot(x,nanmean(y2,1),'Color',c2);
    if doFill==1
        fill([x fliplr(x)],[nanmean(y2,1)+nanstd(y2,[],1)./sqrt(size(y2,1)) fliplr(nanmean(y2,1)-nanstd(y2,[],1)./sqrt(size(y2,1)))],[0.1 0.7 0.5]);
    end
    plot(x,nanmean(y2,1),'Color',c2);
    plot(x,nanmean(y2,1)+nanstd(y2,[],1)./sqrt(size(y2,1)),'Color',c2);
    plot(x,nanmean(y2,1)-nanstd(y2,[],1)./sqrt(size(y2,1)),'Color',c2);
    
    plot(x,nanmean(y1,1),'Color',c1);
end
    
end

function powerInBand=takeTimeSlice(specgram,timeWindow,noTheta)

if length(size(specgram))>2
    powerInBand=nan(size(specgram,3),size(specgram,2)); 
    for i=1:size(specgram,3)
        powerInBand(i,:)=reshape(nanmean(specgram(noTheta.allS.t>=timeWindow(1) & noTheta.allS.t<=timeWindow(2),:,i),1),size(specgram,2),1);
    end
else
    powerInBand=nan(1,size(specgram,2)); 
    powerInBand(1,:)=reshape(nanmean(specgram(noTheta.allS.t>=timeWindow(1) & noTheta.allS.t<=timeWindow(2),:),1),size(specgram,2),1);
end

end

function powerInBand=takeFrequencyBand(specgram,freqBand,noTheta)

if length(size(specgram))>2
    powerInBand=nan(size(specgram,3),size(specgram,1)); 
    for i=1:size(specgram,3)
        powerInBand(i,:)=reshape(nanmean(specgram(:,noTheta.allS.f>=freqBand(1) & noTheta.allS.f<=freqBand(2),i),2),size(specgram,1),1);
    end
else
    powerInBand=nan(1,size(specgram,1)); 
    powerInBand(1,:)=reshape(nanmean(specgram(:,noTheta.allS.f>=freqBand(1) & noTheta.allS.f<=freqBand(2)),2),size(specgram,1),1);
end

end

function combined=addSigMask(data,sigMask,noTheta,p)

tempfill=nanmean(nanmean(data,1),2);
combined=data.*(sigMask(noTheta.allS.t<12,noTheta.allS.f<=50)<p);
combined(combined==0)=tempfill;

end

function sigMask=calculateSignificanceMask(dataset1,dataset2)

% dim 3 is units

sigMask=nan(size(dataset1,1),size(dataset2,2));
for i=1:size(dataset1,1)
    for j=1:size(dataset1,2)
        temp1=dataset1(i,j,:);
        temp2=dataset2(i,j,:);
        if all(isnan(temp1)) || all(isnan(temp2))
            sigMask(i,j)=1;
        else
            sigMask(i,j)=signrank(reshape(temp1,size(temp1,3),1),reshape(temp2,size(temp2,3),1));        
        end
    end
end

end

function output=normAll(input)

for i=1:size(input,3)
    temp=reshape(input(:,:,i),size(input,1),size(input,2));
    normTo=nansum(temp(:,5));
    temp(:,1)=(temp(:,1)./nansum(temp(:,1))).*normTo;
    temp(:,2)=(temp(:,2)./nansum(temp(:,2))).*normTo;
    temp(:,3)=(temp(:,3)./nansum(temp(:,3))).*normTo;
    temp(:,4)=(temp(:,4)./nansum(temp(:,4))).*normTo;
    temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
    temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
    output(:,:,i)=temp;
end

end

function output=whitenAll(input,isCx,isHigh)

if isCx==false
    for i=1:size(input,3)
        temp=reshape(input(:,:,i),size(input,1),size(input,2));
        temp(:,1:3)=temp(:,1:3)./(1.9*0.17.*(0.7225/0.1121));
        temp(:,4)=temp(:,4)./(1.4*0.73.*(0.7225/0.5035));
        temp(:,2)=1.07*temp(:,2);
        temp(:,4)=1.07*temp(:,4);
        temp(:,3)=1.0475*temp(:,3);
        temp(:,1:3)=temp(:,1:3)+0.1;
        temp(:,4)=temp(:,4)+0.1;
        leaveMinAt=nanmin(temp(:,3));
        temp(:,3)=temp(:,3)*1.5-1.5*leaveMinAt-0.1;
        temp(:,2)=temp(:,2)*1.2-0.9*leaveMinAt+0.18;
        temp(:,4)=temp(:,4)*1.2-0.9*leaveMinAt+0.18;
        output(:,:,i)=temp;
    end
%     for i=1:size(input,3)
%         temp=reshape(input(:,:,i),size(input,1),size(input,2));
%         temp(:,1:3)=temp(:,1:3)./(1.9*0.17.*(0.7225/0.1121));
%         temp(:,4)=temp(:,4)./(1.4*0.73.*(0.7225/0.5035));
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,4)=1.07*temp(:,4);
%         temp(:,3)=1.0475*temp(:,3);
%         output(:,:,i)=temp;
%     end
else
    for i=1:size(input,3)
%         temp=reshape(input(:,:,i),size(input,1),size(input,2));
%         temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*0.37);
%         temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.39);
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,1)=1.05*temp(:,1);
%         temp(:,3)=1.06*temp(:,3);
%         temp(:,4)=0.75*temp(:,4);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
%         output(:,:,i)=temp;
         if isHigh==0
            temp=reshape(input(:,:,i),size(input,1),size(input,2));
            temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*1*1.1);
            temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.7*1.1);
            temp(:,2)=0.95*temp(:,2);
            temp(:,1)=1.0*temp(:,1);
            temp(:,3)=1.07*temp(:,3);
            temp(:,4)=0.42*temp(:,4);
            temp(:,5)=0.48*temp(:,5);
            temp(:,6)=0.95*temp(:,6);
            output(:,:,i)=temp;
            
         else
%             temp=reshape(input(:,:,i),size(input,1),size(input,2));
%             temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*1.5);
%             temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.9);
%             temp(:,2)=1.07*temp(:,2);
%             temp(:,1)=1.03*temp(:,1);
%             temp(:,3)=1.06*temp(:,3);
%             temp(:,4)=0.75*temp(:,4);
%             output(:,:,i)=temp;
            
            temp=reshape(input(:,:,i),size(input,1),size(input,2));
            temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*1*1.3);
            temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.7*1.3);
            temp(:,2)=0.95*temp(:,2);
            temp(:,1)=1.0*temp(:,1);
            temp(:,3)=1.07*temp(:,3);
            temp(:,4)=0.42*temp(:,4);
            temp(:,5)=0.43*temp(:,5);
            temp(:,6)=0.95*temp(:,6);
            output(:,:,i)=temp;
         end
    end
end

end

function temp=whitenAndNormSpecgram(temp,isCx,isHigh)

if isCx==false
    temp(:,1:3)=temp(:,1:3)./(1.9*0.17.*(0.7225/0.1121));
    temp(:,4)=temp(:,4)./(1.4*0.73.*(0.7225/0.5035));
    temp(:,2)=1.07*temp(:,2);
    temp(:,4)=1.07*temp(:,4);
    temp(:,3)=1.0475*temp(:,3);
    % temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
    % temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
else
    if isHigh==0
%         temp(:,1:3)=temp(:,1:3)./(3.3*0.17.*(0.7225/0.1121));
%         temp(:,4)=temp(:,4)./(2.0*0.73.*(0.7225/0.5035));
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,4)=1.07*temp(:,4);
%         %     temp(:,3)=1.0475*temp(:,3);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));


%         % BEST FOR NO THETA
%         temp(:,1:3)=temp(:,1:3)./(3.53*0.17.*(0.7225/0.1121));
%         temp(:,4)=temp(:,4)./(2.13*0.73.*(0.7225/0.5035));
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,4)=1.07*temp(:,4);
%         %     temp(:,3)=1.0475*temp(:,3);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));


        temp(:,1:3)=temp(:,1:3)./(4*0.17.*(0.7225/0.1121));
        temp(:,4)=temp(:,4)./(3*0.73.*(0.7225/0.5035));
        temp(:,2)=1.07*temp(:,2);
        temp(:,4)=1.07*temp(:,4);
        %     temp(:,3)=1.0475*temp(:,3);
        temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
        temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
    else
%         temp(:,1:3)=temp(:,1:3)./(3.3*0.17.*(0.7225/0.1121)*1.2);
%         temp(:,4)=temp(:,4)./(2.0*0.73.*(0.7225/0.5035)*0.9);
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,4)=1.03*temp(:,4);
%         temp(:,3)=1.0*temp(:,3);
%         temp(:,4)=0.75*temp(:,4);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));

%         % BEST FOR NO THETA
%         temp(:,1:3)=temp(:,1:3)./(3.53*0.17.*(0.7225/0.1121));
%         temp(:,4)=temp(:,4)./(2.13*0.73.*(0.7225/0.5035));
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,4)=1.07*temp(:,4);
%         %     temp(:,3)=1.0475*temp(:,3);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));

        temp(:,1:3)=temp(:,1:3)./(4*0.17.*(0.7225/0.1121));
        temp(:,4)=temp(:,4)./(3*0.73.*(0.7225/0.5035));
        temp(:,2)=1.07*temp(:,2);
        temp(:,4)=1.07*temp(:,4);
        %     temp(:,3)=1.0475*temp(:,3);
        temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
        temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
    end
end

K=ones([2 2]); 
temp=conv2(temp,K,'same');
temp=interp2(temp,4);

temp=temp(1:end-13,1:end-13);

end

function out=addAll(data1,data2)

tmp=cat(3,data1,data2);
out=nansum(tmp,3);

end

function align_specgram=alignToStim(specgram,specgram_t,placeStimOnsetAt,currStimOnset)

if placeStimOnsetAt==currStimOnset
    align_specgram=specgram;
    return
end

timeDiff=abs(currStimOnset-placeStimOnsetAt);
inds_diff=floor(timeDiff/(specgram_t(2)-specgram_t(1)));

if currStimOnset>placeStimOnsetAt
    align_specgram=[specgram(inds_diff:end,:); nan(inds_diff-1,size(specgram,2))];
elseif currStimOnset<placeStimOnsetAt
    align_specgram=[nan(inds_diff,size(specgram,2)); specgram(1:end-inds_diff,:)];
end

end

function a=loadVar(varargin)

if length(varargin)==3
    varDir=varargin{1};
    varName=varargin{2};
    loadFieldByField=varargin{3};
    inFileName=[];
elseif length(varargin)==4
    varDir=varargin{1};
    varName=varargin{2};
    loadFieldByField=varargin{3};
    inFileName=varargin{4};
else
    error('expected 3 or 4 arguments to loadVar');
end

if loadFieldByField==true
    if isempty(inFileName)
        a.(varName)=loadStructFieldByField([varDir '\' varName]);
        if iscell(a.(varName))
            a=a.(varName);
        end
    else
        a.(inFileName)=loadStructFieldByField([varDir '\' varName]);
        if iscell(a.(inFileName))
            a=a.(inFileName);
        end
    end
else
    a=load([varDir '\' varName '.mat']);
end

end
       