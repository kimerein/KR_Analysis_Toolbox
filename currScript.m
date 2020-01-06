tookFileInds=[1:155];
trialDuration=14.5;
outputDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Ntsr1 opto stim\Mawake397\dLGN T1 LED sines';
ledOffVal=[0 0.05];
% take=[];
putStimAt=4; % in seconds

detectfilenames=expt.files.names(tookFileInds);
spikes=fixSpikesAfterO2Sort(spikes,expt,detectfilenames);
save([outputDir '\' spikes.info.spikesfile '.mat'],'spikes');
temp=spikes.labels(:,1);
useAssigns=temp(spikes.labels(:,2)==2); % 2 means good unit
disp(useAssigns);
psth=measureAllUnitsResponseProperties(spikes,useAssigns,[0 trialDuration]);
save([outputDir '\psth_notAlignedToStim.mat'],'psth');
psth=alignToStimFromExpt(psth,expt,tookFileInds,trialDuration,putStimAt);
save([outputDir '\psth_alignedToStim.mat'],'psth');
s=unique(spikes.stimcond(~isnan(spikes.stimcond)));
disp(s);
clear uses; clear uses_tri;
if isempty(take)
    a=unique(spikes.trials(~isnan(spikes.stimcond))); % use all trials
else
    a=take; % use only these trials
end
for i=1:length(s)
    uses{i}={[s(i)]};
    uses_tri{i}={a};
end
led=unique(psth.unitLED{1});
ledOnVal=led(led>=0.01);
doF1analysis([],[],outputDir,[],psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials);











% for i=1:length(dLGNpsth.psths)
%     temp=dLGNpsth.psths{i};
%     temp=downSampMatrix(temp,2);
%     dLGNpsth.psths{i}=temp;    
% end
% dLGNpsth.t=downSampAv(dLGNpsth.t,2);



% function currScript(psth,noThetaTrials)
% 
% outputDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Simultaneous dLGN and Ntsr1\Mawake401\dLGN\specs';
% 
% V1dir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Simultaneous dLGN and Ntsr1\Mawake401\V1 Ntsr1\specs';
% 
% s=unique(psth.unitStimcond{1});
% clear uses; clear uses_tri;
% a=psth.unitTrials{1};
% for i=1:length(s)
%     uses{i}={[s(i)]};
%     uses_tri{i}={a};
% end
% ledOffVal=[0 0.05];
% led=unique(psth.unitLED{1});
% ledOnVal=led(~ismember(led,ledOffVal));
% ledOnVal=100; % to save time
% 
% % Comment out if already done
% % doF1analysis_saveAllTrialSpecs([],[],outputDir,[],psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials);
% 
% load([outputDir '\noTheta_trialAv_temp_noLED.mat']);
% [lowF1Trials,highF1Trials]=findTrialsWithHighOrLowF1(noTheta_trialAv_temp,psth,uses,uses_tri,noThetaTrials,1);
% load([outputDir '\noTheta_noLED.mat']);
% firststep=find(noThetaTrials'==1 & ismember(psth.unitLED{1},ledOffVal));
% indsIntoSpecs=find(ismember(lowF1Trials,firststep));
% temp=zeros(size(noTheta.allS.S{1}));
% for i=1:length(noTheta.allS.S)
%     temp=temp+noTheta.allS.S{i};
% end
% temp=temp./length(noTheta.allS.S);
% temp=reshape(nanmean(temp(:,:,indsIntoSpecs),3),size(temp,1),size(temp,2));
% figure(); 
% imagesc(noTheta.allS.t,noTheta.allS.f,temp');
% title('No theta low F1');
% indsIntoSpecs=find(ismember(highF1Trials,firststep));
% temp=noTheta.allS.S{1};
% temp=reshape(nanmean(temp(:,:,indsIntoSpecs),3),size(temp,1),size(temp,2));
% figure(); 
% imagesc(noTheta.allS.t,noTheta.allS.f,temp');
% title('No theta high F1');
% 
% % Plot V1
% load([V1dir '\noTheta_noLED.mat']);
% firststep=find(noThetaTrials'==1 & ismember(psth.unitLED{1},ledOffVal));
% temp=zeros(size(noTheta.allS.S{1}));
% for i=1:length(noTheta.allS.S)
%     temp=temp+noTheta.allS.S{i};
% end
% temp=temp./length(noTheta.allS.S);
% temp=reshape(nanmean(temp(:,:,indsIntoSpecs),3),size(temp,1),size(temp,2));
% figure(); 
% imagesc(noTheta.allS.t,noTheta.allS.f,temp');
% title('V1 units -- No theta low F1');
% indsIntoSpecs=find(ismember(highF1Trials,firststep));
% temp=noTheta.allS.S{1};
% temp=reshape(nanmean(temp(:,:,indsIntoSpecs),3),size(temp,1),size(temp,2));
% figure(); 
% imagesc(noTheta.allS.t,noTheta.allS.f,temp');
% title('V1 units -- No theta high F1');
% 
% 
% load([outputDir '\theta_trialAv_temp_noLED.mat']);
% [lowF1Trials,highF1Trials]=findTrialsWithHighOrLowF1(theta_trialAv_temp,psth,uses,uses_tri,noThetaTrials,0);
% load([outputDir '\theta_noLED.mat']);
% firststep=find(noThetaTrials'==0 & ismember(psth.unitLED{1},ledOffVal));
% indsIntoSpecs=find(ismember(lowF1Trials,firststep));
% temp=zeros(size(theta.allS.S{1}));
% for i=1:length(theta.allS.S)
%     temp=temp+theta.allS.S{i};
% end
% temp=temp./length(theta.allS.S);
% temp=reshape(nanmean(temp(:,:,indsIntoSpecs),3),size(temp,1),size(temp,2));
% figure(); 
% imagesc(theta.allS.t,theta.allS.f,temp');
% title('Theta low F1');
% indsIntoSpecs=find(ismember(highF1Trials,firststep));
% temp=theta.allS.S{1};
% temp=reshape(nanmean(temp(:,:,indsIntoSpecs),3),size(temp,1),size(temp,2));
% figure(); 
% imagesc(theta.allS.t,theta.allS.f,temp');
% title('Theta high F1');
% 
% % Plot V1
% load([V1dir '\theta_noLED.mat']);
% indsIntoSpecs=find(ismember(lowF1Trials,firststep));
% temp=zeros(size(theta.allS.S{1}));
% for i=1:length(theta.allS.S)
%     temp=temp+theta.allS.S{i};
% end
% temp=temp./length(theta.allS.S);
% temp=reshape(nanmean(temp(:,:,indsIntoSpecs),3),size(temp,1),size(temp,2));
% figure(); 
% imagesc(theta.allS.t,theta.allS.f,temp');
% title('V1 units -- Theta low F1');
% indsIntoSpecs=find(ismember(highF1Trials,firststep));
% temp=theta.allS.S{1};
% temp=reshape(nanmean(temp(:,:,indsIntoSpecs),3),size(temp,1),size(temp,2));
% figure(); 
% imagesc(theta.allS.t,theta.allS.f,temp');
% title('V1 units -- Theta high F1');


% 
% % for i=1:size(all_forP_noTheta_spont,2)
% %     for j=1:size(all_forP_noTheta_spont,3)
% %         temper1=all_forP_noTheta_spont(:,i,j);
% %         temper2=all_forP_noTheta_visEv(:,i,j);
% %         p_noTheta(i,j)=signrank(temper1(1:end),temper2(1:end));
% %         temper1=all_forP_theta_spont(:,i,j);
% %         temper2=all_forP_theta_visEv(:,i,j);
% %         p_theta(i,j)=signrank(temper1(1:end),temper2(1:end));
% %     end
% % end
% 
% use_noTheta_specs=noTheta_allS(allClassifyAsNtsr1==1);
% use_theta_specs=theta_allS(allClassifyAsNtsr1==1);
% 
% clear allExVals_noTheta
% clear allExVals_theta
% clear allExSpikes_noTheta
% clear allExSpikes_theta
% allExVals_noTheta=nan(1,length(use_noTheta_specs));
% allExVals_theta=nan(1,length(use_noTheta_specs));
% allAmps_noTheta=nan(1,length(use_noTheta_specs));
% allAmps_theta=nan(1,length(use_noTheta_specs));
% allExSpikes_noTheta=nan(1,length(use_noTheta_specs));
% allExSpikes_theta=nan(1,length(use_noTheta_specs));
% use_f=f(f<=50);
% for i=1:length(use_noTheta_specs)
%     tempSpec=use_noTheta_specs{i};
% %     tempSpec=use_noTheta_specs{i}+use_theta_specs{i};
%     forP_noTheta_visEv(i,:)=nanmean(tempSpec(times>1 & times<=3,:),1);
%     forP_noTheta_spont(i,:)=nanmean(tempSpec(times<1,:),1);
% %     noTheta_visR=nanmean(tempSpec(times>1 & times<=3,:),1)-nanmean(tempSpec(times<1,:),1);
%     noTheta_visR=nanmean(tempSpec(times>1.5 & times<=3,:),1);
% %     noTheta_visR=smooth(noTheta_visR,2);
% %     noTheta_visR(1:3)=noTheta_visR(1:3).*0.17.*(0.7225/0.1121);
% %     noTheta_visR(4)=noTheta_visR(4).*0.73.*(0.7225/0.5035);
%     noTheta_visR(1:3)=noTheta_visR(1:3).*1.2*0.17.*(0.7225/0.1121);
%     noTheta_visR(4)=noTheta_visR(4).*1.2*0.73.*(0.7225/0.5035);
%     highAlphaPower=nanmean(noTheta_visR(f>9.5 & f<20));
%     highAlphaAmp=sqrt(abs(highAlphaPower));
% %     highAlphaAmp=sign(highAlphaPower)*highAlphaAmp;
%     lowAlphaPower=nanmean(noTheta_visR(f>9.5 & f<20));
%     lowAlphaAmp=sqrt(abs(lowAlphaPower));
% %     lowAlphaAmp=sign(lowAlphaPower)*lowAlphaAmp;
% %     noTheta_visR(f<4 & f>30)=nanmin(noTheta_visR);
%     noTheta_visR=noTheta_visR(f<=50);
% %     noTheta_visR=noTheta_visR-nanmin(noTheta_visR);
% %     noTheta_visR=noTheta_visR./nanmax(noTheta_visR);
% %     noTheta_expectedVal=nansum((noTheta_visR./nansum(noTheta_visR)).*use_f);
%     [~,ma]=nanmax(noTheta_visR);
%     allExSpikes_noTheta(i)=use_f(ma).*sqrt(nanmean(noTheta_visR(use_f>use_f(ma)-1 & use_f<use_f(ma)+1)));
%     allAmps_noTheta(i)=sqrt(nanmean(noTheta_visR(use_f>use_f(ma)-1 & use_f<use_f(ma)+1)));
%     allExVals_noTheta(i)=use_f(ma);
%     
%     tempSpec=use_theta_specs{i};
%     forP_theta_visEv(i,:)=nanmean(tempSpec(times>1 & times<=3,:),1);
%     forP_theta_spont(i,:)=nanmean(tempSpec(times<1,:),1);
% %     noTheta_visR=nanmean(tempSpec(times>1 & times<=3,:),1)-nanmean(tempSpec(times<1,:),1);
%     noTheta_visR=nanmean(tempSpec(times>1.5 & times<=3,:),1);
% %     noTheta_visR=smooth(noTheta_visR,2);
% %     noTheta_visR(1:3)=noTheta_visR(1:3).*0.17.*(0.7225/0.1121);
% %     noTheta_visR(4)=noTheta_visR(4).*0.73.*(0.7225/0.5035);
%     noTheta_visR(1:3)=noTheta_visR(1:3).*1.2*0.17.*(0.7225/0.1121);
%     noTheta_visR(4)=noTheta_visR(4).*1.2*0.73.*(0.7225/0.5035);
%     highAlphaPower=nanmean(noTheta_visR(f>9.5 & f<20));
%     highAlphaAmp=sqrt(abs(highAlphaPower));
% %     highAlphaAmp=sign(highAlphaPower)*highAlphaAmp;
%     lowAlphaPower=nanmean(noTheta_visR(f>9.5 & f<20));
%     lowAlphaAmp=sqrt(abs(lowAlphaPower));
% %     lowAlphaAmp=sign(lowAlphaPower)*lowAlphaAmp;
% %     noTheta_visR(f<4 & f>30)=nanmin(noTheta_visR);
%     noTheta_visR=noTheta_visR(f<=50);
% %     noTheta_visR=noTheta_visR-nanmin(noTheta_visR);
% %     noTheta_visR=noTheta_visR./nanmax(noTheta_visR);
% %     theta_expectedVal=nansum((noTheta_visR./nansum(noTheta_visR)).*use_f);
%     [~,ma]=nanmax(noTheta_visR);
%     allExSpikes_theta(i)=use_f(ma).*sqrt(nanmean(noTheta_visR(use_f>use_f(ma)-1 & use_f<use_f(ma)+1)));
%     allAmps_theta(i)=sqrt(nanmean(noTheta_visR(use_f>use_f(ma)-1 & use_f<use_f(ma)+1)));
%     allExVals_theta(i)=use_f(ma);
% end
% 
% k=15;
% allFreqs_exSpikes_noTheta(k,:)=allExSpikes_noTheta;
% allFreqs_exSpikes_theta(k,:)=allExSpikes_theta;
% 
% allPeakFreqs_noTheta(k,:)=allExVals_noTheta;
% allPeakFreqs_theta(k,:)=allExVals_theta;
% 
% allFreqs_amps_noTheta(k,:)=allAmps_noTheta;
% allFreqs_amps_theta(k,:)=allAmps_theta;
% 
% % for j=1:length(freqs)
% %     
% % 
% % plotThisFreq=freqs(j);
% % dd{1}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Ntsr1 opto tagging\AR2019-08-15 V1\FF\OneDrive-2019-11-04\T1_noAv' '\Hz' num2str(plotThisFreq)];
% % dd{2}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Ntsr1 opto tagging\AR2019-08-15 V1\FF\OneDrive-2019-11-04\T2_noAv' '\Hz' num2str(plotThisFreq)];
% % dd{3}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Ntsr1 opto tagging\AR2019-08-15 V1\FF\OneDrive-2019-11-04\T3_noAv' '\Hz' num2str(plotThisFreq)];
% % dd{4}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Ntsr1 opto tagging\AR2019-08-15 V1\FF\OneDrive-2019-11-04\T4_noAv' '\Hz' num2str(plotThisFreq)];
% % dd{5}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Ntsr1 opto tagging\AR2019-08-15 V1\FF\OneDrive-2019-11-04\T5_noAv' '\Hz' num2str(plotThisFreq)];
% % dd{6}=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Ntsr1 opto tagging\AR2019-08-15 V1\FF\OneDrive-2019-11-04\T6_noAv' '\Hz' num2str(plotThisFreq)];
% % 
% % sub.(['output_' num2str(plotThisFreq) 'Hz'])=plotF1analysis_simple(dd,trialDuration,times,allLines);
% % close all;
% % i=find(ismember(freqs,plotThisFreq)); 
% % trytemp=sub.(['output_' num2str(plotThisFreq) 'Hz']);
% % noTheta_HFa(i,:)=nanmean(trytemp.noTheta_noLED(:,times>1.5 & times<=2.5),2); theta_HFa(i,:)=nanmean(trytemp.theta_noLED(:,times>1.5 & times<=2.5),2);
% % 
% % end