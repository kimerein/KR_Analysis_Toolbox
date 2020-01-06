function output=plotF1analysis_simple(datadir,trialDuration,times,allLines)

if iscell(datadir)
    all_noTheta_trialAv_noLED=[];
    all_noTheta_trialAv_LED=[];
    all_theta_trialAv_noLED=[];
    all_theta_trialAv_LED=[];
    
    for i=1:length(datadir)
        d=datadir{i};
%         a=load([d '\' 'noTheta_trialAv_noLED']);
%         noTheta_trialAv_noLED=a.noTheta_trialAv;
        a=load([d '\' 'noTheta_noLED']);
        noTheta_trialAv_noLED=a.noTheta;
        all_noTheta_trialAv_noLED=concatStructs(all_noTheta_trialAv_noLED,noTheta_trialAv_noLED);
        
%         a=load([d '\' 'noTheta_trialAv_LED']);
%         noTheta_trialAv_LED=a.noTheta_trialAv;
        a=load([d '\' 'noTheta_LED']);
        noTheta_trialAv_LED=a.noTheta;
        all_noTheta_trialAv_LED=concatStructs(all_noTheta_trialAv_LED,noTheta_trialAv_LED);
        
%         a=load([d '\' 'theta_trialAv_noLED']);
%         theta_trialAv_noLED=a.theta_trialAv;
        a=load([d '\' 'theta_noLED']);
        theta_trialAv_noLED=a.theta;
        all_theta_trialAv_noLED=concatStructs(all_theta_trialAv_noLED,theta_trialAv_noLED);

%         a=load([d '\' 'theta_trialAv_LED']);
%         theta_trialAv_LED=a.theta_trialAv;
        a=load([d '\' 'theta_LED']);
        theta_trialAv_LED=a.theta;
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

% plotWStderr_simple(noTheta_trialAv_noLED.F1amp,theta_trialAv_noLED.F1amp,trialDuration,'k','r',times,allLines);
% [noTheta_noLED,noTheta_LED,isBigEnough]=plotWStderr_simple(noTheta_trialAv_noLED.F1amp,noTheta_trialAv_LED.F1amp,trialDuration,'k','b',times,allLines);
% [theta_noLED,theta_LED,isBigEnough]=plotWStderr_simple(theta_trialAv_noLED.F1amp,theta_trialAv_LED.F1amp,trialDuration,'r','c',times,allLines);
% plotWStderr_simple(noTheta_trialAv_noLED.F1amp,theta_trialAv_LED.F1amp,trialDuration,'k','c',times,allLines);

plotWStderr_simple(noTheta_trialAv_noLED.HFa,theta_trialAv_noLED.HFa,trialDuration,'k','r',times,allLines);
[noTheta_noLED,noTheta_LED,isBigEnough]=plotWStderr_simple(noTheta_trialAv_noLED.HFa,noTheta_trialAv_LED.HFa,trialDuration,'k','b',times,allLines);
[theta_noLED,theta_LED,isBigEnough]=plotWStderr_simple(theta_trialAv_noLED.HFa,theta_trialAv_LED.HFa,trialDuration,'r','c',times,allLines);
plotWStderr_simple(noTheta_trialAv_noLED.HFa,theta_trialAv_LED.HFa,trialDuration,'k','c',times,allLines);

% plotWStderr_simple(noTheta_trialAv_noLED.LFa,theta_trialAv_noLED.LFa,trialDuration,'k','r',times,allLines);
% [noTheta_noLED,noTheta_LED,isBigEnough]=plotWStderr_simple(noTheta_trialAv_noLED.LFa,noTheta_trialAv_LED.LFa,trialDuration,'k','b',times,allLines);
% [theta_noLED,theta_LED,isBigEnough]=plotWStderr_simple(theta_trialAv_noLED.LFa,theta_trialAv_LED.LFa,trialDuration,'r','c',times,allLines);
% plotWStderr_simple(noTheta_trialAv_noLED.LFa,theta_trialAv_LED.LFa,trialDuration,'k','c',times,allLines);

output.noTheta_noLED=noTheta_noLED;
output.noTheta_LED=noTheta_LED;
output.theta_noLED=theta_noLED;
output.theta_LED=theta_LED;

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