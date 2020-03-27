function getDenominatorForFreqResponse(filedir,stimTimeWindow,spontTimeWindow)

doAsAmp=false;

filedir=[filedir '\'];
a=load([filedir 'noTheta_trialAv_temp_noLED.mat']);
noTheta_noLED=a.noTheta_trialAv_temp;
a=load([filedir 'theta_trialAv_temp_noLED.mat']);
theta_noLED=a.theta_trialAv_temp;
a=load([filedir 'noTheta_trialAv_temp_LED.mat']);
noTheta_LED=a.noTheta_trialAv_temp;
a=load([filedir 'theta_trialAv_temp_LED.mat']);
theta_LED=a.theta_trialAv_temp;

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% make 15 x 15 grid
% for each unit
units_noTheta_noLED.p=cell(1,length(noTheta_noLED(1).allS.S));
units_noTheta_LED.p=cell(1,length(noTheta_noLED(1).allS.S));
units_theta_noLED.p=cell(1,length(noTheta_noLED(1).allS.S));
units_theta_LED.p=cell(1,length(noTheta_noLED(1).allS.S));
units_noTheta_noLED.pStim=cell(1,length(noTheta_noLED(1).allS.S));
units_noTheta_LED.pStim=cell(1,length(noTheta_noLED(1).allS.S));
units_theta_noLED.pStim=cell(1,length(noTheta_noLED(1).allS.S));
units_theta_LED.pStim=cell(1,length(noTheta_noLED(1).allS.S));
for j=1:length(units_noTheta_noLED.p)
    currGrid_noTheta_noLED=nan(15,15);
    currGrid_noTheta_LED=nan(15,15);
    currGrid_theta_noLED=nan(15,15);
    currGrid_theta_LED=nan(15,15);
    currGrid_noTheta_noLED_notBaseSub=nan(15,15);
    currGrid_noTheta_LED_notBaseSub=nan(15,15);
    currGrid_theta_noLED_notBaseSub=nan(15,15);
    currGrid_theta_LED_notBaseSub=nan(15,15);
    for i=1:length(freqs)
        curr_trialAv=noTheta_noLED(i);
        S=curr_trialAv.allS.S{j};
        t=curr_trialAv.allS.t;
        f=curr_trialAv.allS.f;
        for k=1:length(freqs)
            [resp,nonBaseSub_resp]=getResponse(S,t,f,stimTimeWindow,spontTimeWindow,freqs(k),doAsAmp);
            currGrid_noTheta_noLED(i,k)=resp;
            currGrid_noTheta_noLED_notBaseSub(i,k)=nonBaseSub_resp;
        end
        curr_trialAv=theta_noLED(i);
        S=curr_trialAv.allS.S{j};
        t=curr_trialAv.allS.t;
        f=curr_trialAv.allS.f;
        for k=1:length(freqs)
            [resp,nonBaseSub_resp]=getResponse(S,t,f,stimTimeWindow,spontTimeWindow,freqs(k),doAsAmp);
            currGrid_theta_noLED(i,k)=resp;
            currGrid_noTheta_LED_notBaseSub(i,k)=nonBaseSub_resp;
        end
        
        curr_trialAv=noTheta_LED(i);
        S=curr_trialAv.allS.S{j};
        t=curr_trialAv.allS.t;
        f=curr_trialAv.allS.f;
        for k=1:length(freqs)
            [resp,nonBaseSub_resp]=getResponse(S,t,f,stimTimeWindow,spontTimeWindow,freqs(k),doAsAmp);
            currGrid_noTheta_LED(i,k)=resp;
            currGrid_theta_noLED_notBaseSub(i,k)=nonBaseSub_resp;
        end
        
        curr_trialAv=theta_LED(i);
        S=curr_trialAv.allS.S{j};
        t=curr_trialAv.allS.t;
        f=curr_trialAv.allS.f;
        for k=1:length(freqs)
            [resp,nonBaseSub_resp]=getResponse(S,t,f,stimTimeWindow,spontTimeWindow,freqs(k),doAsAmp);
            currGrid_theta_LED(i,k)=resp;
            currGrid_theta_LED_notBaseSub(i,k)=nonBaseSub_resp;
        end
    end
    units_noTheta_noLED.p{j}=currGrid_noTheta_noLED;
    units_noTheta_LED.p{j}=currGrid_noTheta_LED;
    units_theta_noLED.p{j}=currGrid_theta_noLED;
    units_theta_LED.p{j}=currGrid_theta_LED;
    
    units_noTheta_noLED.pStim{j}=currGrid_noTheta_noLED_notBaseSub;
    units_noTheta_LED.pStim{j}=currGrid_noTheta_LED_notBaseSub;
    units_theta_noLED.pStim{j}=currGrid_theta_noLED_notBaseSub;
    units_theta_LED.pStim{j}=currGrid_theta_LED_notBaseSub;
end

% Get diagonal and integral for each cell
diags_units_noTheta_noLED=nan(length(noTheta_noLED(1).allS.S),15);
diags_units_noTheta_LED=nan(length(noTheta_noLED(1).allS.S),15);
diags_units_theta_noLED=nan(length(noTheta_noLED(1).allS.S),15);
diags_units_theta_LED=nan(length(noTheta_noLED(1).allS.S),15);
ints_noTheta_noLED=nan(1,length(noTheta_noLED(1).allS.S));
ints_noTheta_LED=nan(1,length(noTheta_noLED(1).allS.S));
ints_theta_noLED=nan(1,length(noTheta_noLED(1).allS.S));
ints_theta_LED=nan(1,length(noTheta_noLED(1).allS.S));
for j=1:length(units_noTheta_noLED.p)
    for i=1:15
        temp=units_noTheta_noLED.p{j};
        diags_units_noTheta_noLED(j,i)=temp(i,i);
        if i==1
            temp=units_noTheta_noLED.p{j};
            ints_noTheta_noLED(j)=nansum(nansum(temp));
        end
        temp=units_noTheta_LED.p{j};
        diags_units_noTheta_LED(j,i)=temp(i,i);
        if i==1
            temp=units_noTheta_LED.p{j};
            ints_noTheta_LED(j)=nansum(nansum(temp));
        end
        temp=units_theta_noLED.p{j};
        diags_units_theta_noLED(j,i)=temp(i,i);
        if i==1
            temp=units_theta_noLED.p{j};
            ints_theta_noLED(j)=nansum(nansum(temp));
        end
        temp=units_theta_LED.p{j};
        diags_units_theta_LED(j,i)=temp(i,i);
        if i==1
            temp=units_theta_LED.p{j};
            ints_theta_LED(j)=nansum(nansum(temp));
        end
    end
end

save([filedir 'diags_units_noTheta_noLED.mat'],'diags_units_noTheta_noLED');
save([filedir 'diags_units_noTheta_LED.mat'],'diags_units_noTheta_LED');
save([filedir 'diags_units_theta_noLED.mat'],'diags_units_theta_noLED');
save([filedir 'diags_units_theta_LED.mat'],'diags_units_theta_LED');
save([filedir 'ints_noTheta_noLED.mat'],'ints_noTheta_noLED');
save([filedir 'ints_noTheta_LED.mat'],'ints_noTheta_LED');
save([filedir 'ints_theta_noLED.mat'],'ints_theta_noLED');
save([filedir 'ints_theta_LED.mat'],'ints_theta_LED');

end

function [resp,nonBaseSub_resp]=getResponse(S,t,f,stimTimeWindow,spontTimeWindow,freq,doAsAmp)

% Get F1 power stim minus spont
stim=nanmean(nanmean(S(t>stimTimeWindow(1) & t<=stimTimeWindow(2),f>=freq-0.75 & f<=freq+0.75),1),2);
spont=nanmean(nanmean(S(t>spontTimeWindow(1) & t<=spontTimeWindow(2),f>=freq-0.75 & f<=freq+0.75),1),2);
resp=stim-spont;
nonBaseSub_resp=stim;
    
% Change to amplitude?
if doAsAmp==true
    resp=sqrt(resp);
    nonBaseSub_resp=sqrt(nonBaseSub_resp);
end

end

