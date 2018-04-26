%getUnitsFRsAcrossTrials

fileInd=3:70;
useAssigns=[4 31 39 43];
T='1';
dataDir='W:\Analysis Computer\Persistent Units\Mawake64\';
useLED=[3 5];
useStimcond=1:32;

sspikes=filtspikes(spikes,0,'fileInd',fileInd);
subSpikes=filtspikes(sspikes,0,'led',useLED,'stimcond',useStimcond);

% pers=cell(1,length(useAssigns));
% before=cell(1,length(useAssigns));
for i=1:length(useAssigns)
%       [m,s,n]=calcMAndSForUnit(filtspikes(sspikes,0,'assigns',useAssigns(i)),[0.933 1.04],unique(subSpikes.trials));
%     [m,s,n]=calcMAndSForUnit(filtspikes(sspikes,0,'assigns',useAssigns(i)),[0.893 1],unique(subSpikes.trials));
%     [m,s,n]=calcMAndSForUnit(filtspikes(sspikes,0,'assigns',useAssigns(i)),[1.303 1.41],unique(subSpikes.trials));
    [m,s,n]=calcMAndSForUnit(filtspikes(sspikes,0,'assigns',useAssigns(i)),[1.343 1.45],unique(subSpikes.trials));
%     pers{i}=n;
%       [mb,sb,nb]=calcMAndSForUnit(filtspikes(sspikes,0,'assigns',useAssigns(i)),[0 0.107],unique(subSpikes.trials));
%     [mb,sb,nb]=calcMAndSForUnit(filtspikes(sspikes,0,'assigns',useAssigns(i)),[0 0.107],unique(subSpikes.trials));
%     [mb,sb,nb]=calcMAndSForUnit(filtspikes(sspikes,0,'assigns',useAssigns(i)),[1.193 1.30],unique(subSpikes.trials));
    [mb,sb,nb]=calcMAndSForUnit(filtspikes(sspikes,0,'assigns',useAssigns(i)),[1.233 1.34],unique(subSpikes.trials));
%     before{i}=nb;
    save([dataDir 'T_' T '_' num2str(useAssigns(i)) 'pers35.mat'],'n');
    save([dataDir 'T_' T '_' num2str(useAssigns(i)) 'before35.mat'],'nb');
end
