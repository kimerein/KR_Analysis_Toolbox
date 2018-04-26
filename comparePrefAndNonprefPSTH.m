function [varargout]=comparePrefAndNonprefPSTH(spikes,useTheseAssigns,useTheseStimcond,bin,ledOnCond)

onsetResponseWindow=[1.05 1.3];
useThisLEDwindow=[1.26 1.32];
% onsetResponseWindow=[3.05 3.3];
% useThisLEDwindow=[3.36 3.42];

prefR=zeros(length(useTheseAssigns),1);
prefStimcond=zeros(length(useTheseAssigns),1);
nonprefR=zeros(length(useTheseAssigns),1);
nonprefStimcond=zeros(length(useTheseAssigns),1);
sigStim=zeros(length(useTheseAssigns),1);
for i=1:length(useTheseAssigns)
    stimcondR=zeros(length(useTheseStimcond),1);
    for j=1:length(useTheseStimcond)
        [stimcondR(j),temp,d{j}]=calcMeanAndStdDuringWindow(filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',useTheseStimcond(j)),onsetResponseWindow);
    end
    [prefR(i),ind]=max(stimcondR);
    prefStimcond(i)=useTheseStimcond(ind);
    [nonprefR(i),ind]=min(stimcondR);
    nonprefStimcond(i)=useTheseStimcond(ind);
    sigStim(i)=mattest(d{prefStimcond(i)}',d{nonprefStimcond(i)}');
end
disp('Done with calculating preferred and non-preferred responses');

sigDiff=zeros(length(useTheseAssigns),1);
prefStimLED=zeros(length(useTheseAssigns),1);
nonprefStimLED=zeros(length(useTheseAssigns),1);
for i=1:length(useTheseAssigns)
    disp(i);
    if i==1
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',prefStimcond(i));
        [n1,centers1,edges1,xpoints1,ypoints1,stds1]=psth_wStdev_valuesOnly(filtSpikes(useSpikes,0,'led',ledOnCond),bin);
        prefPSTH_xpoints=xpoints1;
        prefPSTH_ypoints=ypoints1;
        [prefStimLED(i),s,prefR_LED]=calcMeanAndStdDuringWindow(useSpikes,useThisLEDwindow);
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',nonprefStimcond(i));
        [n1,centers1,edges1,xpoints1,ypoints1,stds1]=psth_wStdev_valuesOnly(filtSpikes(useSpikes,0,'led',ledOnCond),bin);
        nonprefPSTH_xpoints=xpoints1;
        nonprefPSTH_ypoints=ypoints1;
        [nonprefStimLED(i),s,nonprefR_LED]=calcMeanAndStdDuringWindow(useSpikes,useThisLEDwindow);
%         disp(size(prefR_LED));
%         disp(size(nonprefR_LED));
        sigDiff(i)=mattest(prefR_LED',nonprefR_LED');
    else
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',prefStimcond(i));
        [n1,centers1,edges1,xpoints1,ypoints1,stds1]=psth_wStdev_valuesOnly(filtSpikes(useSpikes,0,'led',ledOnCond),bin);
        prefPSTH_ypoints=prefPSTH_ypoints+ypoints1;
        [prefStimLED(i),s,prefR_LED]=calcMeanAndStdDuringWindow(useSpikes,useThisLEDwindow);
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',nonprefStimcond(i));
        [n1,centers1,edges1,xpoints1,ypoints1,stds1]=psth_wStdev_valuesOnly(filtSpikes(useSpikes,0,'led',ledOnCond),bin);
        nonprefPSTH_ypoints=nonprefPSTH_ypoints+ypoints1;
        [nonprefStimLED(i),s,nonprefR_LED]=calcMeanAndStdDuringWindow(useSpikes,useThisLEDwindow);
%         disp(size(prefR_LED));
%         disp(size(nonprefR_LED));
        sigDiff(i)=mattest(prefR_LED',nonprefR_LED');
    end
end
figure(); 
plot(prefPSTH_xpoints,prefPSTH_ypoints./length(useTheseAssigns),'Color','blue');
hold on;
plot(nonprefPSTH_xpoints,nonprefPSTH_ypoints./length(useTheseAssigns),'Color',[0.2 0.2 0.2]);
[h,p]=ttest(prefStimLED,nonprefStimLED,0.05,'both');

varargout{1}=h;
varargout{2}=p;
varargout{3}=prefStimLED;
varargout{4}=nonprefStimLED;
varargout{5}=sigDiff;
varargout{6}=prefStimcond;
varargout{7}=nonprefStimcond;
varargout{8}=prefR;
varargout{9}=nonprefR;
varargout{10}=sigStim;
varargout{11}=prefPSTH_xpoints;
varargout{12}=prefPSTH_ypoints./length(useTheseAssigns);
varargout{13}=nonprefPSTH_ypoints./length(useTheseAssigns);

    