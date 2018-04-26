function [prefAv,nonprefAv,allrandAv,shuf_prefAv,shuf_nonprefAv,pref_trialSynch,nonpref_trialSynch]=synchCalc_wrapper(spikes,useAssigns,orderCells,refFreq)

trialDuration=14.5;
% trialDuration=8;
synchWindow=15;
s=1;
% s=1:4;
% s=[2 4];
usel=[0];

allrandAv=[];

prefAv=zeros(length(s),floor(trialDuration/(synchWindow/1000)));
disp('Doing pref av');
for i=1:length(s)
    disp(i);
    [prefAv(i,:),pref_trialSynch{i}]=plotAlphaRasters(spikes,useAssigns,orderCells,[],[],0,refFreq,s(i));
    close all
end
shuf_prefAv=zeros(length(s),floor(trialDuration/(synchWindow/1000)));
disp('Doing trial-shuffled pref av');
permd_spikes=spikes;
for i=1:length(s)
    disp(i);
    currTrials=unique(spikes.sweeps.trials(ismember(spikes.sweeps.led,usel) & ismember(spikes.sweeps.stimcond,s(i))));
    % Shuffle trials for each unit separately
    a=unique(spikes.assigns);
    for k=1:length(a)
        randMapping=currTrials(randperm(length(currTrials)));
        switched=zeros(1,length(permd_spikes.trials));
        for j=1:length(currTrials)
            takeTrials=permd_spikes.assigns==a(k) & permd_spikes.trials==currTrials(j) & switched==0;
            permd_spikes.trials(takeTrials)=randMapping(j);
            switched(takeTrials)=1;
        end
    end
    shuf_prefAv(i,:)=plotAlphaRasters(permd_spikes,useAssigns,orderCells,[],[],0,refFreq,i);
end
    

nonprefAv=zeros(length(s),floor(trialDuration/(synchWindow/1000)));
disp('Doing nonpref av');
for i=1:length(s)
    disp(i);
    [nonprefAv(i,:),nonpref_trialSynch{i}]=plotAlphaRasters(spikes,useAssigns,-orderCells,[],[],0,refFreq,s(i));
    close all
end
shuf_nonprefAv=zeros(length(s),floor(trialDuration/(synchWindow/1000)));
disp('Doing trial-shuffled nonpref av');
for i=1:length(s)
    disp(i);
    shuf_nonprefAv(i,:)=plotAlphaRasters(permd_spikes,useAssigns,-orderCells,[],[],0,refFreq,s(i));
end

% % Bootstrap random unit sets
% n=1;
% allrandAv=zeros(length(s)*n,floor(trialDuration/(synchWindow/1000)));
% disp('Doing bootstrap rand of n trials');
% for j=1:n
%     disp(j);
%     randAv=zeros(length(s),floor(trialDuration/(synchWindow/1000)));
%     for i=1:length(s)
%         r=nan(size(orderCells));
%         for k=1:size(orderCells,2)
%             r(:,k)=orderCells(randperm(size(orderCells,1)),k);
%         end
%         randAv(i,:)=plotAlphaRasters(spikes,useAssigns,r,[],[],0,refFreq,i);
%         close all
%     end
%     allrandAv((j-1)*length(s)+1:j*length(s),:)=randAv;
% end

figure(); 
plot(linspace(0,trialDuration,floor(trialDuration/(synchWindow/1000))),nanmean(prefAv,1),'Color','r');
hold on;
plot(linspace(0,trialDuration,floor(trialDuration/(synchWindow/1000))),nanmean(prefAv,1)+nanstd(prefAv,[],1)./sqrt(size(prefAv,1)),'Color','r');
plot(linspace(0,trialDuration,floor(trialDuration/(synchWindow/1000))),nanmean(prefAv,1)-nanstd(prefAv,[],1)./sqrt(size(prefAv,1)),'Color','r');

plot(linspace(0,trialDuration,floor(trialDuration/(synchWindow/1000))),nanmean(nonprefAv,1),'Color','b');
plot(linspace(0,trialDuration,floor(trialDuration/(synchWindow/1000))),nanmean(nonprefAv,1)+nanstd(nonprefAv,[],1)./sqrt(size(nonprefAv,1)),'Color','b');
plot(linspace(0,trialDuration,floor(trialDuration/(synchWindow/1000))),nanmean(nonprefAv,1)-nanstd(nonprefAv,[],1)./sqrt(size(nonprefAv,1)),'Color','b');


% plot(linspace(0,trialDuration,floor(trialDuration/(synchWindow/1000))),nanmean(allrandAv,1),'Color','k');
% plot(linspace(0,trialDuration,floor(trialDuration/(synchWindow/1000))),nanmean(allrandAv,1)+nanstd(allrandAv,[],1)./sqrt(size(allrandAv,1)),'Color','k');
% plot(linspace(0,trialDuration,floor(trialDuration/(synchWindow/1000))),nanmean(allrandAv,1)-nanstd(allrandAv,[],1)./sqrt(size(allrandAv,1)),'Color','k');