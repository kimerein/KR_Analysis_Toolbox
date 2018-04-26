

% temp=filtspikes(s234,0,'fileInd',90:163); 
% 
% [~,~,~,x,y,n]=psth_wStd_trialByTrial(temp,1,0,3.5,length(unique(temp.sweeps.trials)),unique(temp.sweeps.trials));

% clear newt; 
% t=n;
% for i=1:size(t,1)
% newt(i,:)=downSampAv(t(i,:),1);
% end

newt=t;
ttemp{1}=newt; 

s=getAlphaSpikes(expt,filtspikes(s234,0,'fileInd',90:163),ttemp,[]);