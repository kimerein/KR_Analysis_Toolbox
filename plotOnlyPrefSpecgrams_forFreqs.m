function plotOnlyPrefSpecgrams_forFreqs(d,outdir,norm)

% F1range=[2.5 3.5];
% stimWindow=[4 6.5];
% stimWindow=[4.4 6.5];
% spontWindow=[0 4];
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

stimWindow=[1 3];
% stimWindow=[2 3];
% spontWindow=[0.5 1]; % for LOW FREQ
spontWindow=[0 0.7]; % for HIGH FREQ

a=load([d '\' 'noTheta_trialAv_temp_noLED']);
noTheta_trialAv_noLED=a.noTheta_trialAv_temp;

a=load([d '\' 'noTheta_trialAv_temp_LED']);
noTheta_trialAv_LED=a.noTheta_trialAv_temp;

a=load([d '\' 'theta_trialAv_temp_noLED']);
theta_trialAv_noLED=a.theta_trialAv_temp;

a=load([d '\' 'theta_trialAv_temp_LED']);
theta_trialAv_LED=a.theta_trialAv_temp;

t=theta_trialAv_noLED(1).allS.t;
f=theta_trialAv_noLED(1).allS.f;
prefSpecgrams_noTheta_noLED=cell(1,length(theta_trialAv_noLED(1).allS.S));
prefSpecgrams_noTheta_LED=cell(1,length(theta_trialAv_noLED(1).allS.S));
prefSpecgrams_theta_noLED=cell(1,length(theta_trialAv_noLED(1).allS.S));
prefSpecgrams_theta_LED=cell(1,length(theta_trialAv_noLED(1).allS.S));
nonprefSpecgrams_noTheta_noLED=cell(1,length(theta_trialAv_noLED(1).allS.S));
nonprefSpecgrams_noTheta_LED=cell(1,length(theta_trialAv_noLED(1).allS.S));
nonprefSpecgrams_theta_noLED=cell(1,length(theta_trialAv_noLED(1).allS.S));
nonprefSpecgrams_theta_LED=cell(1,length(theta_trialAv_noLED(1).allS.S));

prefSpecgrams_noTheta_noLED_freqs=cell(1,length(theta_trialAv_noLED(1).allS.S));
prefSpecgrams_noTheta_LED_freqs=cell(1,length(theta_trialAv_noLED(1).allS.S));
prefSpecgrams_theta_noLED_freqs=cell(1,length(theta_trialAv_noLED(1).allS.S));
prefSpecgrams_theta_LED_freqs=cell(1,length(theta_trialAv_noLED(1).allS.S));
nonprefSpecgrams_noTheta_noLED_freqs=cell(1,length(theta_trialAv_noLED(1).allS.S));
nonprefSpecgrams_noTheta_LED_freqs=cell(1,length(theta_trialAv_noLED(1).allS.S));
nonprefSpecgrams_theta_noLED_freqs=cell(1,length(theta_trialAv_noLED(1).allS.S));
nonprefSpecgrams_theta_LED_freqs=cell(1,length(theta_trialAv_noLED(1).allS.S));
for i=1:length(theta_trialAv_noLED(1).allS.S)
    F1power_theta_noLED=nan(1,length(theta_trialAv_noLED));
    F1power_theta_LED=nan(1,length(theta_trialAv_noLED));
    F1power_noTheta_LED=nan(1,length(theta_trialAv_noLED));
    F1power_noTheta_noLED=nan(1,length(theta_trialAv_noLED));
%     for j=2:4
%     for j=11:12
    for j=1:length(theta_trialAv_noLED)
        F1range(1)=freqs(j)-0.5;
        F1range(2)=freqs(j)+0.5;
        if j==1
            F1range(1)=0.5;
            F1range(2)=1.6;
        end
        currS=theta_trialAv_noLED(j).allS.S{i};
        meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
        F1=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),f>=F1range(1) & f<=F1range(2)),1),2)-nanmean(nanmean(meancurrS(t>=spontWindow(1) & t<=spontWindow(2),f>=F1range(1) & f<=F1range(2)),1),2);
%         F1=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),f>=F1range(1) & f<=F1range(2)),1),2);
        if norm==0
            F1power_theta_noLED(j)=F1;
        else
            allpower=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),:),1),2);
            F1power_theta_noLED(j)=F1./allpower;
        end
        currS=noTheta_trialAv_LED(j).allS.S{i};
        meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
        F1=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),f>=F1range(1) & f<=F1range(2)),1),2)-nanmean(nanmean(meancurrS(t>=spontWindow(1) & t<=spontWindow(2),f>=F1range(1) & f<=F1range(2)),1),2);
%         F1=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),f>=F1range(1) & f<=F1range(2)),1),2);
        if norm==0
            F1power_noTheta_LED(j)=F1;
        else
            allpower=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),:),1),2);
            F1power_noTheta_LED(j)=F1./allpower;
        end
        currS=noTheta_trialAv_noLED(j).allS.S{i};
        meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
        F1=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),f>=F1range(1) & f<=F1range(2)),1),2)-nanmean(nanmean(meancurrS(t>=spontWindow(1) & t<=spontWindow(2),f>=F1range(1) & f<=F1range(2)),1),2);
%         F1=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),f>=F1range(1) & f<=F1range(2)),1),2);
        if norm==0
            F1power_noTheta_noLED(j)=F1;
        else
            allpower=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),:),1),2);
            F1power_noTheta_noLED(j)=F1./allpower;
        end
        currS=theta_trialAv_LED(j).allS.S{i};
        meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
        F1=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),f>=F1range(1) & f<=F1range(2)),1),2)-nanmean(nanmean(meancurrS(t>=spontWindow(1) & t<=spontWindow(2),f>=F1range(1) & f<=F1range(2)),1),2);
%         F1=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),f>=F1range(1) & f<=F1range(2)),1),2);
        if norm==0
            F1power_theta_LED(j)=F1;
        else
            allpower=nanmean(nanmean(meancurrS(t>=stimWindow(1) & t<=stimWindow(2),:),1),2);
            F1power_theta_LED(j)=F1./allpower;
        end
    end
%     meanF1=nanmean([F1power_theta_noLED; F1power_noTheta_LED; F1power_theta_LED],1);
%     meanF1=nanmean([F1power_theta_noLED; F1power_theta_LED],1);
      meanF1=nanmean([F1power_theta_noLED],1);
%     meanF1=nanmean([F1power_noTheta_LED],1);
%   meanF1=nanmean([F1power_noTheta_noLED],1);
%     meanF1=nanmean([F1power_theta_LED],1);
%     meanF1=nanmean([F1power_noTheta_LED],1);
    [s,si]=sort(meanF1,'descend');
    si=si(~isnan(s));
    prefStim=si(1);
    nonprefStim=si(end);
    
    currS=noTheta_trialAv_noLED(prefStim).allS.S{i};
    meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
    prefSpecgrams_noTheta_noLED{i}=meancurrS;
    prefSpecgrams_noTheta_noLED_freqs{i}=freqs(prefStim);
    
    currS=noTheta_trialAv_LED(prefStim).allS.S{i};
    meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
    prefSpecgrams_noTheta_LED{i}=meancurrS;
    prefSpecgrams_noTheta_LED_freqs{i}=freqs(prefStim);
    
    currS=theta_trialAv_noLED(prefStim).allS.S{i};
    meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
    prefSpecgrams_theta_noLED{i}=meancurrS;
    prefSpecgrams_theta_noLED_freqs{i}=freqs(prefStim);
    
    currS=theta_trialAv_LED(prefStim).allS.S{i};
    meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
    prefSpecgrams_theta_LED{i}=meancurrS;
    prefSpecgrams_theta_LED_freqs{i}=freqs(prefStim);
    
    
    
    currS=noTheta_trialAv_noLED(nonprefStim).allS.S{i};
    meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
    nonprefSpecgrams_noTheta_noLED{i}=meancurrS;
    nonprefSpecgrams_noTheta_noLED_freqs{i}=freqs(nonprefStim);
    
    currS=noTheta_trialAv_LED(nonprefStim).allS.S{i};
    meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
    nonprefSpecgrams_noTheta_LED{i}=meancurrS;
    nonprefSpecgrams_noTheta_LED_freqs{i}=freqs(nonprefStim);
    
    currS=theta_trialAv_noLED(nonprefStim).allS.S{i};
    meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
    nonprefSpecgrams_theta_noLED{i}=meancurrS;
    nonprefSpecgrams_theta_noLED_freqs{i}=freqs(nonprefStim);
    
    currS=theta_trialAv_LED(nonprefStim).allS.S{i};
    meancurrS=reshape(nanmean(currS,3),size(currS,1),size(currS,2));
    nonprefSpecgrams_theta_LED{i}=meancurrS;
    nonprefSpecgrams_theta_LED_freqs{i}=freqs(nonprefStim);
end
       

save([outdir '\' 'prefSpecgrams_noTheta_noLED.mat'],'prefSpecgrams_noTheta_noLED');
save([outdir '\' 'prefSpecgrams_noTheta_LED.mat'],'prefSpecgrams_noTheta_LED');
save([outdir '\' 'prefSpecgrams_theta_noLED.mat'],'prefSpecgrams_theta_noLED');
save([outdir '\' 'prefSpecgrams_theta_LED.mat'],'prefSpecgrams_theta_LED');

save([outdir '\' 'prefSpecgrams_noTheta_noLED_freqs.mat'],'prefSpecgrams_noTheta_noLED_freqs');
save([outdir '\' 'prefSpecgrams_noTheta_LED_freqs.mat'],'prefSpecgrams_noTheta_LED_freqs');
save([outdir '\' 'prefSpecgrams_theta_noLED_freqs.mat'],'prefSpecgrams_theta_noLED_freqs');
save([outdir '\' 'prefSpecgrams_theta_LED_freqs.mat'],'prefSpecgrams_theta_LED_freqs');

save([outdir '\' 'nonprefSpecgrams_noTheta_noLED.mat'],'nonprefSpecgrams_noTheta_noLED');
save([outdir '\' 'nonprefSpecgrams_noTheta_LED.mat'],'nonprefSpecgrams_noTheta_LED');
save([outdir '\' 'nonprefSpecgrams_theta_noLED.mat'],'nonprefSpecgrams_theta_noLED');
save([outdir '\' 'nonprefSpecgrams_theta_LED.mat'],'nonprefSpecgrams_theta_LED');

save([outdir '\' 'nonprefSpecgrams_noTheta_noLED_freqs.mat'],'nonprefSpecgrams_noTheta_noLED_freqs');
save([outdir '\' 'nonprefSpecgrams_noTheta_LED_freqs.mat'],'nonprefSpecgrams_noTheta_LED_freqs');
save([outdir '\' 'nonprefSpecgrams_theta_noLED_freqs.mat'],'nonprefSpecgrams_theta_noLED_freqs');
save([outdir '\' 'nonprefSpecgrams_theta_LED_freqs.mat'],'nonprefSpecgrams_theta_LED_freqs');

save([outdir '\' 't.mat'],'t');
save([outdir '\' 'f.mat'],'f');