function currForStats=combinePrefSpecgrams_forFreqs(datadir,trialDuration)

% F1range=[2.5 3.5];
% stimWindow=[4 6.5];
stimWindow=[1 3];
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

if iscell(datadir)
    all_nonprefSpecgrams_noTheta_LED={};
    all_nonprefSpecgrams_noTheta_noLED={};
    all_nonprefSpecgrams_theta_LED={};
    all_nonprefSpecgrams_theta_noLED={};
    all_nonpref_freqs={};
    
    all_prefSpecgrams_noTheta_LED={};
    all_prefSpecgrams_noTheta_noLED={};
    all_prefSpecgrams_theta_LED={};
    all_prefSpecgrams_theta_noLED={};
    all_pref_freqs={};
    for i=1:length(datadir)
        d=datadir{i};
        a=load([d '\' 't']);
        t=a.t;
        
        d=datadir{i};
        a=load([d '\' 'f']);
        f=a.f;
        
        a=load([d '\' 'nonprefSpecgrams_noTheta_LED']);
        nonprefSpecgrams_noTheta_LED=a.nonprefSpecgrams_noTheta_LED;
        all_nonprefSpecgrams_noTheta_LED=[all_nonprefSpecgrams_noTheta_LED nonprefSpecgrams_noTheta_LED];
        
        a=load([d '\' 'nonprefSpecgrams_noTheta_LED_freqs']);
        nonpref_freq=a.nonprefSpecgrams_noTheta_LED_freqs;
        all_nonpref_freqs=[all_nonpref_freqs nonpref_freq];
        
        a=load([d '\' 'nonprefSpecgrams_noTheta_noLED']);
        nonprefSpecgrams_noTheta_noLED=a.nonprefSpecgrams_noTheta_noLED;
        all_nonprefSpecgrams_noTheta_noLED=[all_nonprefSpecgrams_noTheta_noLED nonprefSpecgrams_noTheta_noLED];
        
        a=load([d '\' 'nonprefSpecgrams_theta_LED']);
        nonprefSpecgrams_theta_LED=a.nonprefSpecgrams_theta_LED;
        all_nonprefSpecgrams_theta_LED=[all_nonprefSpecgrams_theta_LED nonprefSpecgrams_theta_LED];

        a=load([d '\' 'nonprefSpecgrams_theta_noLED']);
        nonprefSpecgrams_theta_noLED=a.nonprefSpecgrams_theta_noLED;
        all_nonprefSpecgrams_theta_noLED=[all_nonprefSpecgrams_theta_noLED nonprefSpecgrams_theta_noLED];
        
        
        a=load([d '\' 'prefSpecgrams_noTheta_LED']);
        prefSpecgrams_noTheta_LED=a.prefSpecgrams_noTheta_LED;
        all_prefSpecgrams_noTheta_LED=[all_prefSpecgrams_noTheta_LED prefSpecgrams_noTheta_LED];
        
        a=load([d '\' 'prefSpecgrams_noTheta_LED_freqs']);
        pref_freq=a.prefSpecgrams_noTheta_LED_freqs;
        all_pref_freqs=[all_pref_freqs pref_freq];
        
        a=load([d '\' 'prefSpecgrams_noTheta_noLED']);
        prefSpecgrams_noTheta_noLED=a.prefSpecgrams_noTheta_noLED;
        all_prefSpecgrams_noTheta_noLED=[all_prefSpecgrams_noTheta_noLED prefSpecgrams_noTheta_noLED];
        
        a=load([d '\' 'prefSpecgrams_theta_LED']);
        prefSpecgrams_theta_LED=a.prefSpecgrams_theta_LED;
        all_prefSpecgrams_theta_LED=[all_prefSpecgrams_theta_LED prefSpecgrams_theta_LED];

        a=load([d '\' 'prefSpecgrams_theta_noLED']);
        prefSpecgrams_theta_noLED=a.prefSpecgrams_theta_noLED;
        all_prefSpecgrams_theta_noLED=[all_prefSpecgrams_theta_noLED prefSpecgrams_theta_noLED];
    end
    nonprefSpecgrams_noTheta_LED=all_nonprefSpecgrams_noTheta_LED;
    nonprefSpecgrams_noTheta_noLED=all_nonprefSpecgrams_noTheta_noLED;
    nonprefSpecgrams_theta_LED=all_nonprefSpecgrams_theta_LED;
    nonprefSpecgrams_theta_noLED=all_nonprefSpecgrams_theta_noLED;
    nonpref_freq=all_nonpref_freqs;
    
    prefSpecgrams_noTheta_LED=all_prefSpecgrams_noTheta_LED;
    prefSpecgrams_noTheta_noLED=all_prefSpecgrams_noTheta_noLED;
    prefSpecgrams_theta_LED=all_prefSpecgrams_theta_LED;
    prefSpecgrams_theta_noLED=all_prefSpecgrams_theta_noLED;
    pref_freq=all_pref_freqs;
else
    d=datadir;
    a=load([d '\' 't']);
    t=a.t;
    
    a=load([d '\' 'f']);
    f=a.f;
    
    a=load([d '\' 'nonprefSpecgrams_noTheta_LED']);
    nonprefSpecgrams_noTheta_LED=a.nonprefSpecgrams_noTheta_LED;
    
    a=load([d '\' 'nonprefSpecgrams_noTheta_noLED']);
    nonprefSpecgrams_noTheta_noLED=a.nonprefSpecgrams_noTheta_noLED;
    
    a=load([d '\' 'nonprefSpecgrams_theta_LED']);
    nonprefSpecgrams_theta_LED=a.nonprefSpecgrams_theta_LED;
    
    a=load([d '\' 'nonprefSpecgrams_theta_noLED']);
    nonprefSpecgrams_theta_noLED=a.nonprefSpecgrams_theta_noLED;
    
    a=load([d '\' 'nonprefSpecgrams_noTheta_LED_freqs']);
    nonpref_freq=a.nonprefSpecgrams_noTheta_LED_freqs;
    
    
    a=load([d '\' 'prefSpecgrams_noTheta_LED']);
    prefSpecgrams_noTheta_LED=a.prefSpecgrams_noTheta_LED;
    
    a=load([d '\' 'prefSpecgrams_noTheta_noLED']);
    prefSpecgrams_noTheta_noLED=a.prefSpecgrams_noTheta_noLED;
    
    a=load([d '\' 'prefSpecgrams_theta_LED']);
    prefSpecgrams_theta_LED=a.prefSpecgrams_theta_LED;
    
    a=load([d '\' 'prefSpecgrams_theta_noLED']);
    prefSpecgrams_theta_noLED=a.prefSpecgrams_theta_noLED;

    a=load([d '\' 'prefSpecgrams_noTheta_LED_freqs']);
    pref_freq=a.prefSpecgrams_noTheta_LED_freqs;
end

% Get F1 and average specgrams
nonprefSpecgrams_noTheta_LED_F1=nan(length(nonprefSpecgrams_noTheta_LED),length(t));
nonprefSpecgrams_noTheta_noLED_F1=nan(length(nonprefSpecgrams_noTheta_LED),length(t));
nonprefSpecgrams_theta_LED_F1=nan(length(nonprefSpecgrams_noTheta_LED),length(t));
nonprefSpecgrams_theta_noLED_F1=nan(length(nonprefSpecgrams_noTheta_LED),length(t));

prefSpecgrams_noTheta_LED_F1=nan(length(nonprefSpecgrams_noTheta_LED),length(t));
prefSpecgrams_noTheta_noLED_F1=nan(length(nonprefSpecgrams_noTheta_LED),length(t));
prefSpecgrams_theta_LED_F1=nan(length(nonprefSpecgrams_noTheta_LED),length(t));
prefSpecgrams_theta_noLED_F1=nan(length(nonprefSpecgrams_noTheta_LED),length(t));

sumSpec.nonprefSpecgrams_noTheta_LED=zeros(size(nonprefSpecgrams_noTheta_LED{1}));
sumSpec.nonprefSpecgrams_noTheta_noLED=zeros(size(nonprefSpecgrams_noTheta_LED{1}));
sumSpec.nonprefSpecgrams_theta_LED=zeros(size(nonprefSpecgrams_noTheta_LED{1}));
sumSpec.nonprefSpecgrams_theta_noLED=zeros(size(nonprefSpecgrams_noTheta_LED{1}));

sumSpec.prefSpecgrams_noTheta_LED=zeros(size(nonprefSpecgrams_noTheta_LED{1}));
sumSpec.prefSpecgrams_noTheta_noLED=zeros(size(nonprefSpecgrams_noTheta_LED{1}));
sumSpec.prefSpecgrams_theta_LED=zeros(size(nonprefSpecgrams_noTheta_LED{1}));
sumSpec.prefSpecgrams_theta_noLED=zeros(size(nonprefSpecgrams_noTheta_LED{1}));
for i=1:length(nonprefSpecgrams_noTheta_LED)
    currS=nonprefSpecgrams_noTheta_LED{i};
    currfreq=nonpref_freq{i};
    if currfreq==1
        F1range(1)=0.5;
        F1range(2)=2;
    else
        F1range(1)=currfreq-0.5;
        F1range(2)=currfreq+0.5;
    end
    nonprefSpecgrams_noTheta_LED_F1(i,:)=nanmean(currS(:,f>=F1range(1) & f<=F1range(2)),2)';
    if ~isnan(currS)
        sumSpec.nonprefSpecgrams_noTheta_LED=sumSpec.nonprefSpecgrams_noTheta_LED+currS;
    end
    
    currS=nonprefSpecgrams_noTheta_noLED{i};
    nonprefSpecgrams_noTheta_noLED_F1(i,:)=nanmean(currS(:,f>=F1range(1) & f<=F1range(2)),2)';
    if ~isnan(currS)
        sumSpec.nonprefSpecgrams_noTheta_noLED=sumSpec.nonprefSpecgrams_noTheta_noLED+currS;
    end
    
    currS=nonprefSpecgrams_theta_LED{i};
    nonprefSpecgrams_theta_LED_F1(i,:)=nanmean(currS(:,f>=F1range(1) & f<=F1range(2)),2)';
    if ~isnan(currS)
        sumSpec.nonprefSpecgrams_theta_LED=sumSpec.nonprefSpecgrams_theta_LED+currS;
    end
    
    currS=nonprefSpecgrams_theta_noLED{i};
    nonprefSpecgrams_theta_noLED_F1(i,:)=nanmean(currS(:,f>=F1range(1) & f<=F1range(2)),2)';
    if ~isnan(currS)
        sumSpec.nonprefSpecgrams_theta_noLED=sumSpec.nonprefSpecgrams_theta_noLED+currS;
    end
    
    
    
    currS=prefSpecgrams_noTheta_LED{i};
    currfreq=pref_freq{i};
    if currfreq==1
        F1range(1)=0.5;
        F1range(2)=1.6;
    else
        F1range(1)=currfreq-0.5;
        F1range(2)=currfreq+0.5;
    end
    prefSpecgrams_noTheta_LED_F1(i,:)=nanmean(currS(:,f>=F1range(1) & f<=F1range(2)),2)';
    if ~isnan(currS)
        sumSpec.prefSpecgrams_noTheta_LED=sumSpec.prefSpecgrams_noTheta_LED+currS;
    end
    
    currS=prefSpecgrams_noTheta_noLED{i};
    prefSpecgrams_noTheta_noLED_F1(i,:)=nanmean(currS(:,f>=F1range(1) & f<=F1range(2)),2)';
    if ~isnan(currS)
        sumSpec.prefSpecgrams_noTheta_noLED=sumSpec.prefSpecgrams_noTheta_noLED+currS;
    end
    
    currS=prefSpecgrams_theta_LED{i};
    prefSpecgrams_theta_LED_F1(i,:)=nanmean(currS(:,f>=F1range(1) & f<=F1range(2)),2)';
    if ~isnan(currS)
        sumSpec.prefSpecgrams_theta_LED=sumSpec.prefSpecgrams_theta_LED+currS;
    end
    
    currS=prefSpecgrams_theta_noLED{i};
    prefSpecgrams_theta_noLED_F1(i,:)=nanmean(currS(:,f>=F1range(1) & f<=F1range(2)),2)';
    if ~isnan(currS)
        sumSpec.prefSpecgrams_theta_noLED=sumSpec.prefSpecgrams_theta_noLED+currS;
    end
end

% Plot F1
[nonprefSpecgrams_noTheta_noLED_F1,nonprefSpecgrams_noTheta_LED_F1]=plotWStderr(nonprefSpecgrams_noTheta_noLED_F1,nonprefSpecgrams_noTheta_LED_F1,trialDuration,'b','r');
% thresh=0;
% isBigEnough=nanmean(nonprefSpecgrams_noTheta_noLED_F1(:,t>stimWindow(1) & t<=stimWindow(2)),2)>thresh;
% [nonprefSpecgrams_theta_noLED_F1,nonprefSpecgrams_theta_LED_F1]=plotWStderr(nonprefSpecgrams_theta_noLED_F1(isBigEnough==1),nonprefSpecgrams_theta_LED_F1(isBigEnough==1),trialDuration,'b','c');
% [prefSpecgrams_noTheta_noLED_F1,prefSpecgrams_noTheta_LED_F1]=plotWStderr(prefSpecgrams_noTheta_noLED_F1(isBigEnough==1),prefSpecgrams_noTheta_LED_F1(isBigEnough==1),trialDuration,'b','r');
% [prefSpecgrams_theta_noLED_F1,prefSpecgrams_theta_LED_F1]=plotWStderr(prefSpecgrams_theta_noLED_F1(isBigEnough==1),prefSpecgrams_theta_LED_F1(isBigEnough==1),trialDuration,'b','c');
[nonprefSpecgrams_theta_noLED_F1,nonprefSpecgrams_theta_LED_F1]=plotWStderr(nonprefSpecgrams_theta_noLED_F1,nonprefSpecgrams_theta_LED_F1,trialDuration,'b','c');
[prefSpecgrams_noTheta_noLED_F1,prefSpecgrams_noTheta_LED_F1]=plotWStderr(prefSpecgrams_noTheta_noLED_F1,prefSpecgrams_noTheta_LED_F1,trialDuration,'b','r');
[prefSpecgrams_theta_noLED_F1,prefSpecgrams_theta_LED_F1]=plotWStderr(prefSpecgrams_theta_noLED_F1,prefSpecgrams_theta_LED_F1,trialDuration,'b','c');




% Get statistics on ratio of pref to nonpref
perc1=25;
perc2=75;

% To choose only visually responsive units
% currNonpref_LED=nanmean(nonprefSpecgrams_noTheta_LED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
% currPref_LED=nanmean(prefSpecgrams_noTheta_LED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currNonpref_LED=nanmean(nonprefSpecgrams_theta_noLED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currPref_LED=nanmean(prefSpecgrams_theta_noLED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);

figure(); 
currNonpref=nanmean(nonprefSpecgrams_noTheta_noLED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currPref=nanmean(prefSpecgrams_noTheta_noLED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currPref=currPref(currPref_LED>prctile(currPref_LED,75));
currNonpref=currNonpref(currPref_LED>prctile(currPref_LED,75));
% curr=currPref./currNonpref;
curr=(currPref-currNonpref)./(currPref+currNonpref);
curr(isinf(curr))=nan;
scatter(1,nanmedian(curr));
% scatter(1,nanmean(curr));
hold on;
curr_forStats_noTheta_noLED=curr;
line([1 1],[prctile(curr,perc1) prctile(curr,perc2)]);
% line([1 1],[nanmean(curr)-nanstd(curr,[],1)./sqrt(length(curr)) nanmean(curr)+nanstd(curr,[],1)./sqrt(length(curr))]);

currNonpref=nanmean(nonprefSpecgrams_noTheta_LED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currPref=nanmean(prefSpecgrams_noTheta_LED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currPref=currPref(currPref_LED>prctile(currPref_LED,75));
currNonpref=currNonpref(currPref_LED>prctile(currPref_LED,75));
% curr=currPref./currNonpref;
curr=(currPref-currNonpref)./(currPref+currNonpref);
curr(isinf(curr))=nan;
scatter(2,nanmedian(curr));
% scatter(2,nanmean(curr));
hold on;
curr_forStats_noTheta_LED=curr;
line([2 2],[prctile(curr,perc1) prctile(curr,perc2)]);
% line([2 2],[nanmean(curr)-nanstd(curr,[],1)./sqrt(length(curr)) nanmean(curr)+nanstd(curr,[],1)./sqrt(length(curr))]);

currNonpref=nanmean(nonprefSpecgrams_theta_noLED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currPref=nanmean(prefSpecgrams_theta_noLED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currPref=currPref(currPref_LED>prctile(currPref_LED,75));
currNonpref=currNonpref(currPref_LED>prctile(currPref_LED,75));
% curr=currPref./currNonpref;
curr=(currPref-currNonpref)./(currPref+currNonpref);
curr(isinf(curr))=nan;
scatter(3,nanmedian(curr));
% scatter(3,nanmean(curr));
hold on;
curr_forStats_theta_noLED=curr;
line([3 3],[prctile(curr,perc1) prctile(curr,perc2)]);
% line([3 3],[nanmean(curr)-nanstd(curr,[],1)./sqrt(length(curr)) nanmean(curr)+nanstd(curr,[],1)./sqrt(length(curr))]);

currNonpref=nanmean(nonprefSpecgrams_theta_LED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currPref=nanmean(prefSpecgrams_theta_LED_F1(:,t>=stimWindow(1) & t<=stimWindow(2)),2);
currPref=currPref(currPref_LED>prctile(currPref_LED,75));
currNonpref=currNonpref(currPref_LED>prctile(currPref_LED,75));
% curr=currPref./currNonpref;
curr=(currPref-currNonpref)./(currPref+currNonpref);
curr(isinf(curr))=nan;
scatter(4,nanmedian(curr));
% scatter(4,nanmean(curr));
hold on;
curr_forStats_theta_LED=curr;
line([4 4],[prctile(curr,perc1) prctile(curr,perc2)]);
% line([4 4],[nanmean(curr)-nanstd(curr,[],1)./sqrt(length(curr)) nanmean(curr)+nanstd(curr,[],1)./sqrt(length(curr))]);

currForStats.noTheta_noLED=curr_forStats_noTheta_noLED;
currForStats.noTheta_LED=curr_forStats_noTheta_LED;
currForStats.theta_noLED=curr_forStats_theta_noLED;
currForStats.theta_LED=curr_forStats_theta_LED;

% figure(); 
% boxplot(curr_forStats_noTheta_noLED);
% figure(); 
% boxplot(curr_forStats_noTheta_LED);
% figure(); 
% boxplot(curr_forStats_theta_noLED);
% figure(); 
% boxplot(curr_forStats_theta_LED);

p=signrank(curr_forStats_noTheta_noLED,curr_forStats_noTheta_LED);
disp(p);
p=signrank(curr_forStats_theta_noLED,curr_forStats_theta_LED);
disp(p);

% figure();
% [n,x]=hist(curr_forStats_noTheta_noLED,300);
% plot(x,n,'Color','k');
% hold on;
% [n,x]=hist(curr_forStats_noTheta_LED,300);
% plot(x,n,'Color','c');
% 
% figure();
% [n,x]=hist(curr_forStats_theta_noLED,300);
% plot(x,n,'Color','k');
% hold on;
% [n,x]=hist(curr_forStats_theta_LED,300);
% plot(x,n,'Color','c');

return

% Plot average spectrograms
figure(); 
imagesc(t,f(f<=50),sumSpec.nonprefSpecgrams_noTheta_noLED(:,f<=50)');
title('nonprefSpecgrams_noTheta_noLED');

figure(); 
imagesc(t,f(f<=50),sumSpec.nonprefSpecgrams_noTheta_LED(:,f<=50)');
title('nonprefSpecgrams_noTheta_LED');

figure(); 
imagesc(t,f(f<=50),sumSpec.nonprefSpecgrams_theta_noLED(:,f<=50)');
title('nonprefSpecgrams_theta_noLED');

figure(); 
imagesc(t,f(f<=50),sumSpec.nonprefSpecgrams_theta_LED(:,f<=50)');
title('nonprefSpecgrams_theta_LED');





figure(); 
imagesc(t,f(f<=50),sumSpec.prefSpecgrams_noTheta_noLED(:,f<=50)');
title('prefSpecgrams_noTheta_noLED');

figure(); 
imagesc(t,f(f<=50),sumSpec.prefSpecgrams_noTheta_LED(:,f<=50)');
title('prefSpecgrams_noTheta_LED');

figure(); 
imagesc(t,f(f<=50),sumSpec.prefSpecgrams_theta_noLED(:,f<=50)');
title('prefSpecgrams_theta_noLED');

figure(); 
imagesc(t,f(f<=50),sumSpec.prefSpecgrams_theta_LED(:,f<=50)');
title('prefSpecgrams_theta_LED');



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
