function arrangeReversingGratingsByFreq(filedir)

a=load([filedir 'noTheta_trialAv_temp_noLED.mat']);
noTheta_noLED=a.noTheta_trialAv_temp;
a=load([filedir 'theta_trialAv_temp_noLED.mat']);
theta_noLED=a.theta_trialAv_temp;
a=load([filedir 'noTheta_trialAv_temp_LED.mat']);
noTheta_LED=a.noTheta_trialAv_temp;
a=load([filedir 'theta_trialAv_temp_LED.mat']);
theta_LED=a.theta_trialAv_temp;

stimconds=[1 2 3 4];
freqs=[6 12 6 12]; % in Hz
unique_freqs=unique(freqs);
for i=1:length(unique_freqs)
    currs=stimconds(ismember(freqs,unique_freqs(i)));
    for j=1:length(currs)
        if j==1
            noTheta_trialAv_noLED=noTheta_noLED(i);
            theta_trialAv_noLED=theta_noLED(i);
            noTheta_trialAv_LED=noTheta_LED(i);
            theta_trialAv_LED=theta_LED(i);
        else
            noTheta_trialAv_noLED=addFreqResults(noTheta_trialAv_noLED,noTheta_noLED(i));
            theta_trialAv_noLED=addFreqResults(theta_trialAv_noLED,theta_noLED(i));
            noTheta_trialAv_LED=addFreqResults(noTheta_trialAv_LED,noTheta_LED(i));
            theta_trialAv_LED=addFreqResults(theta_trialAv_LED,theta_LED(i));
        end    
    end
    noTheta_trialAv_noLED=divideFreqResults(noTheta_trialAv_noLED,length(currs));
    theta_trialAv_noLED=divideFreqResults(theta_trialAv_noLED,length(currs));
    noTheta_trialAv_LED=divideFreqResults(noTheta_trialAv_LED,length(currs));
    theta_trialAv_LED=divideFreqResults(theta_trialAv_LED,length(currs));
    currdir=[filedir '\Hz' num2str(unique_freqs(i))];
    mkdir(currdir);
    noTheta_trialAv=noTheta_trialAv_noLED;
    theta_trialAv=theta_trialAv_noLED;
    save([currdir '\noTheta_trialAv_noLED.mat'],'noTheta_trialAv');
    save([currdir '\theta_trialAv_noLED.mat'],'theta_trialAv');
    noTheta_trialAv=noTheta_trialAv_LED;
    theta_trialAv=theta_trialAv_LED;
    save([currdir '\noTheta_trialAv_LED.mat'],'noTheta_trialAv');
    save([currdir '\theta_trialAv_LED.mat'],'theta_trialAv');
end

end

function summed=addFreqResults(existing,new)

summed.HFa=existing.HFa+new.HFa;
summed.LFa=existing.LFa+new.LFa;
summed.F1amp=existing.F1amp+new.F1amp;
summed.allpower=existing.allpower+new.allpower;

summed.allS.t=existing.allS.t;
summed.allS.f=existing.allS.f;
for i=1:length(existing.allS.S)
    summed.allS.S{i}=existing.allS.S{i}+new.allS.S{i};
end

end

function existing=divideFreqResults(existing,n)

existing.HFa=existing.HFa/n;
existing.LFa=existing.LFa/n;
existing.F1amp=existing.F1amp/n;
existing.allpower=existing.allpower/n;

for i=1:length(existing.allS.S)
    existing.allS.S{i}=existing.allS.S{i}/n;
end

end
