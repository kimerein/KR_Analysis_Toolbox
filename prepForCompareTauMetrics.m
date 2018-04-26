function [useX,allled,allcon]=prepForCompareTauMetrics(datas)

useSubtraction=0;

indBefore=20;
indAfter=200;

% Pull out tau sections of data
allcon=zeros(length(datas),indBefore+indAfter);
allled=zeros(length(datas),indBefore+indAfter);
ncon=zeros(length(datas),1);
nled=zeros(length(datas),1);
useX=-indBefore+1:1:indAfter;
useX=useX/1000;
for i=1:length(datas)
    ledOnset=datas(i).stimulusOn(1)+datas(i).stimulusDuration/1000;
    firstInd=find(datas(i).xpoints>=ledOnset,1);
    allcon(i,:)=datas(i).ypoints1(firstInd-indBefore+1:firstInd+indAfter);
    allled(i,:)=datas(i).ypoints2(firstInd-indBefore+1:firstInd+indAfter);
    if useSubtraction==1
        allled(i,:)=mean(allcon(i,:)).*ones(size(allcon(i,:)))-(allcon(i,:)-allled(i,:));
    end  
    ncon(i)=datas(i).numtrials1;
    nled(i)=datas(i).numtrials2;
end