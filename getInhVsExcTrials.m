function getInhVsExcTrials(dLGNpsth,V1psth,noThetaTrials)

ledWindow=[3.2 4];
spontWindow=[0 3.2];

if isempty(noThetaTrials)
    noThetaTrials=ones(1,length(dLGNpsth.unitLED{1}));
end

l=dLGNpsth.unitLED{1};
diffWithLED=nan(1,length(l));

summedRate=zeros(size(dLGNpsth.psths{1}));
for i=1:length(dLGNpsth.psths)
    summedRate=summedRate+dLGNpsth.psths{i};
end

spontRate=nanmean(summedRate(ismember(l,5.05) & noThetaTrials==1,dLGNpsth.t>=spontWindow(1) & dLGNpsth.t<=spontWindow(2)),2);
ledRate=nanmean(summedRate(ismember(l,5.05) & noThetaTrials==1,dLGNpsth.t>=ledWindow(1) & dLGNpsth.t<=ledWindow(2)),2);
diffWithLED=spontRate-ledRate;

[n,xout]=hist(diffWithLED,20);
figure();
plot(xout,n);

ds=5;
takeTrials=find(ismember(l,5.05) & noThetaTrials==1);
figure();
plot(downSampAv(dLGNpsth.t,ds),downSampAv(nanmean(summedRate(takeTrials(diffWithLED>0),:),1),ds),'Color','k');

figure();
plot(downSampAv(dLGNpsth.t,ds),downSampAv(nanmean(summedRate(takeTrials(diffWithLED<=0),:),1),ds),'Color','r');


% Get V1 specgrams based on effect of LED on dLGN activity
movingwin=[1 0.05];
params.tapers=[5 6];
params.Fs=1/(V1psth.t(2)-V1psth.t(1));
params.fpass=[1 50];
params.trialave=0;
allS.t=[];
allS.f=[];
allS.S=cell(length(V1psth.psths),1);
for i=1:length(V1psth.psths)
    disp(i);
    p=V1psth.psths{i};
    [S,t,f]=mtspecgrampb(p',movingwin,params);
    allS.S{i}=S;
    if i==1
        allS.t=t;
        allS.f=f;
    end
end

sumLEDdec=zeros(size(allS.S{1},1),size(allS.S{1},2));
sumLEDinc=zeros(size(allS.S{1},1),size(allS.S{1},2));
for i=1:length(allS.S)
    currS=allS.S{i};
    take_currS=nanmean(currS(:,:,takeTrials(diffWithLED>0)),3);
    use_currS=reshape(take_currS,size(currS,1),size(currS,2));
    sumLEDdec=sumLEDdec+use_currS;
end
for i=1:length(allS.S)
    currS=allS.S{i};
    take_currS=nanmean(currS(:,:,takeTrials(diffWithLED<=0)),3);
    use_currS=reshape(take_currS,size(currS,1),size(currS,2));
    sumLEDinc=sumLEDinc+use_currS;
end

figure();
imagesc(t,f,[sumLEDdec' sumLEDinc']);