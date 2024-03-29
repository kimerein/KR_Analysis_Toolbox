function three_PCAs(dLGNpsth,noThetaTrials,s,l)

% N x TC
% average across all trial types but not across neurons
disp(['using ' num2str(length(dLGNpsth.psths)) ' neurons']);
trialaverage_byneuron=nan(length(dLGNpsth.psths),size(dLGNpsth.psths{1},2));
for i=1:length(dLGNpsth.psths)
    temp=dLGNpsth.psths{i};
    temp=temp';
    trialaverage_byneuron(i,:)=nanmean(temp,2);
end
N_x_TC=trialaverage_byneuron';
[coeff,score,latent,tsquared,explained,mu]=pca(N_x_TC);
n=6;
figure(); 
plot(score(:,1:n));
title('trial-averaged neuron response PC');

% NT x C
% average across neurons and trials but not across trial types
if isempty(s)
    s=dLGNpsth.unitStimcond{1};
end
if isempty(l)
    l=dLGNpsth.unitLED{1};
end
disp('using these stimconds');
disp(unique(s));
disp('using these ledconds');
disp(unique(l));
disp('using these brain states');
disp(unique(noThetaTrials));
[unique_rows,~,unique_conds]=unique([s' l' noThetaTrials],'rows');
conds_per_timepoint=repmat(unique_conds,1,size(dLGNpsth.psths{1},2));
conds_per_timepoint=conds_per_timepoint';
added_temp=zeros(size(dLGNpsth.psths{1}));
for i=1:length(dLGNpsth.psths)
    temp=dLGNpsth.psths{i};
    added_temp=added_temp+temp;
end
added_temp=added_temp./length(dLGNpsth.psths);
averaged_neurons=added_temp';
u=unique(conds_per_timepoint(1,:));
NT_x_C=nan(size(dLGNpsth.psths{1},2),length(u));
for i=1:length(u)
    NT_x_C(:,i)=nanmean(averaged_neurons(:,conds_per_timepoint(1,:)==u(i)),2);
end
[coeff,score,latent,tsquared,explained,mu]=pca(NT_x_C);
% for example, PC_n selects these conditions
n=1;
figure(); 
plot(score(:,n));
title('neuron-averaged condition PC');
disp(['top ' num2str(nansum(coeff(:,n)>prctile(coeff(:,n),90))) ' conditions for PC1']);
disp(unique_rows(coeff(:,n)>prctile(coeff(:,n),90),:));

% T x NC
% average population neural activity vector for each trial
neuronsprofile_withintrial=nan(length(dLGNpsth.psths),size(dLGNpsth.psths{1},1));
for i=1:length(dLGNpsth.psths)
    temp=dLGNpsth.psths{i};
    neuronsprofile_withintrial(i,:)=nanmean(temp,2);
end
T_x_NC=neuronsprofile_withintrial';
[coeff,score,latent,tsquared,explained,mu]=pca(T_x_NC);
n=1;
figure(); 
plot(score(:,n));
% hold on;
% plot(noThetaTrials*10,'Color','k');
temp=corrcoef(score(:,n),s);
disp(['PC' num2str(n) ' corrcoef with stim cond: ' num2str(temp(1,2))]);
temp=corrcoef(score(:,n),l);
disp(['PC' num2str(n) ' corrcoef with led cond: ' num2str(temp(1,2))]);
temp=corrcoef(score(:,n),double(noThetaTrials'));
disp(['PC' num2str(n) ' corrcoef with brain state cond: ' num2str(temp(1,2))]);


end