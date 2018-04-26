function findAlphaUnits(data)

alphaRange=[12 20];

t=data.t;
f=data.f;

alphaPower=zeros(length(data.low.S),1);

for i=1:length(data.low.S)
    curr1=data.low.S{i};
    curr2=data.high.S{i};
    meancurr=(curr1+curr2)./2;
    alphaPower(i)=nanmean(nanmean(meancurr(:,f>=alphaRange(1) & f<=alphaRange(2)),2),1)./nanmean(nanmean(meancurr,2),1);
end

[n,xout]=hist(alphaPower,100);
figure(); 
plot(xout,n);
title('Histogram of Cell Alpha Power');

alphaThresh=1.05;
sumS_low=zeros(size(data.low.S{1}));
sumS_high=zeros(size(data.low.S{1}));
for i=1:length(data.low.S)
    if alphaPower(i)>=alphaThresh
        sumS_low=sumS_low+data.low.S{i};
        sumS_high=sumS_high+data.high.S{i};
    end
end

figure();
imagesc(t,f,[sumS_low' sumS_high']);
    