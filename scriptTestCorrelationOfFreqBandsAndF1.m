tryBands={[1.5 19.5],[1 4],[4 6],[5.5 6.5],[4 8],[6 8],[8 10],[10 12],[11.5 12.5],[12 14],[14 16],[16 18],[9 16],[16 20],[20 40],[40 60]};
    
r=zeros(1,length(tryBands));
p=zeros(1,length(tryBands));

for i=1:length(tryBands)
    
currBand=tryBands{i};
    
% t=trialByTrial_pref{1};
t=pref_trialSynch{1};
params.tapers=[3 10];
params.Fs=1/(14.5/1451);
params.fpass=[1 60];
params.trialave=0;
[S,f]=mtspectrumpb(downSampMatrix(t(:,linspace(0,14.5,1451)>=4 & linspace(0,14.5,1451)<=6.5),1)',params);
[Sbefore,fbefore]=mtspectrumpb(downSampMatrix(t(:,linspace(0,14.5,1451)>=3 & linspace(0,14.5,1451)<=4),1)',params);
x=nanmean(S(f>=2.8 & f<=3.2,:),1);
y=nanmean(Sbefore(fbefore>=currBand(1) & fbefore<=currBand(2),:),1)./nanmean(Sbefore,1);

if currBand(1)==11.5
    figure();
    hist(y,20);
    disp(sum(y>2.5)/length(y));
end

% t=trialByTrial_pref{2};
% params.tapers=[3 10];
% params.Fs=1/(14.5/1451);
% params.fpass=[1 60];
% params.trialave=0;
% [S,f]=mtspectrumpb(downSampMatrix(t(:,linspace(0,14.5,1451)>=4 & linspace(0,14.5,1451)<=6.5),1)',params);
% [Sbefore,fbefore]=mtspectrumpb(downSampMatrix(t(:,linspace(0,14.5,1451)>=3 & linspace(0,14.5,1451)<=5),1)',params);
% x2=nanmean(S(f>=2.8 & f<=3.2,:),1);
% y2=nanmean(Sbefore(fbefore>=currBand(1) & fbefore<=currBand(2),:),1)./nanmean(Sbefore,1);
% 
% 
% t=trialByTrial_pref{3};
% params.tapers=[3 10];
% params.Fs=1/(14.5/1451);
% params.fpass=[1 60];
% params.trialave=0;
% [S,f]=mtspectrumpb(downSampMatrix(t(:,linspace(0,14.5,1451)>=4 & linspace(0,14.5,1451)<=6.5),1)',params);
% [Sbefore,fbefore]=mtspectrumpb(downSampMatrix(t(:,linspace(0,14.5,1451)>=3 & linspace(0,14.5,1451)<=5),1)',params);
% x3=nanmean(S(f>=2.8 & f<=3.2,:),1);
% y3=nanmean(Sbefore(fbefore>=currBand(1) & fbefore<=currBand(2),:),1)./nanmean(Sbefore,1);



% [currr,currp]=corrcoef([x x2 x3],[y y2 y3]);
[currr,currp]=corrcoef([x],[y]);
r(i)=currr(1,2);
p(i)=currp(1,2);
currBandAv(i)=nanmean(currBand);
end

figure(); 
plot(currBandAv,r);

figure(); 
plot(currBandAv,p);
