function plotF1analysis_simple(datadir,trialDuration,times,allLines)

if iscell(datadir)
    all_noTheta_trialAv_noLED=[];
    all_noTheta_trialAv_LED=[];
    all_theta_trialAv_noLED=[];
    all_theta_trialAv_LED=[];
    for i=1:length(datadir)
        d=datadir{i};
        a=load([d '\' 'noTheta_trialAv_noLED']);
        noTheta_trialAv_noLED=a.noTheta_trialAv;
        all_noTheta_trialAv_noLED=concatStructs(all_noTheta_trialAv_noLED,noTheta_trialAv_noLED);
        
        a=load([d '\' 'noTheta_trialAv_LED']);
        noTheta_trialAv_LED=a.noTheta_trialAv;
        all_noTheta_trialAv_LED=concatStructs(all_noTheta_trialAv_LED,noTheta_trialAv_LED);
        
        a=load([d '\' 'theta_trialAv_noLED']);
        theta_trialAv_noLED=a.theta_trialAv;
        all_theta_trialAv_noLED=concatStructs(all_theta_trialAv_noLED,theta_trialAv_noLED);

        a=load([d '\' 'theta_trialAv_LED']);
        theta_trialAv_LED=a.theta_trialAv;
        all_theta_trialAv_LED=concatStructs(all_theta_trialAv_LED,theta_trialAv_LED);
    end
    noTheta_trialAv_noLED=all_noTheta_trialAv_noLED;
    noTheta_trialAv_LED=all_noTheta_trialAv_LED;
    theta_trialAv_noLED=all_theta_trialAv_noLED;
    theta_trialAv_LED=all_theta_trialAv_LED;
elseif ~isempty(datadir)
    a=load([datadir '\' 'noTheta_trialAv_noLED']);
    noTheta_trialAv_noLED=a.noTheta_trialAv;
    
    a=load([datadir '\' 'noTheta_trialAv_LED']);
    noTheta_trialAv_LED=a.noTheta_trialAv;
    
    a=load([datadir '\' 'theta_trialAv_noLED']);
    theta_trialAv_noLED=a.theta_trialAv;
    
    a=load([datadir '\' 'theta_trialAv_LED']);
    theta_trialAv_LED=a.theta_trialAv;
end

plotWStderr_simple(noTheta_trialAv_noLED.F1amp,theta_trialAv_noLED.F1amp,trialDuration,'k','r',times,allLines);
plotWStderr_simple(noTheta_trialAv_noLED.F1amp,noTheta_trialAv_LED.F1amp,trialDuration,'k','b',times,allLines);
plotWStderr_simple(theta_trialAv_noLED.F1amp,theta_trialAv_LED.F1amp,trialDuration,'r','c',times,allLines);

end

function newStruct=concatStructs(struct1,struct2)

if isempty(struct1)
    newStruct=struct2;
    return
elseif isempty(struct2)
    newStruct=struct1;
    return
end
fie=fieldnames(struct1);
for i=1:length(fie)
    f=fie{i};
    if isfield(struct2,f)
        newStruct.(f)=[struct1.(f); struct2.(f)];
    end
end

end