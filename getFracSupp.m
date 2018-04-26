function getFracSupp(datas)



indBefore=20;
indAfter=1000;

% Pull out tau sections of data
allcon=zeros(length(datas),indBefore+indAfter);
allled=zeros(length(datas),indBefore+indAfter);
ncon=zeros(length(datas),1);
nled=zeros(length(datas),1);
useX=-indBefore+1:1:indAfter;
for i=1:length(datas)
    ledOnset=datas(i).stimulusOn(1)+datas(i).stimulusDuration/1000;
    firstInd=find(datas(i).xpoints>=ledOnset,1);
    
    allcon(i,:)=datas(i).ypoints1(firstInd-indBefore+1:firstInd+indAfter);
    allled(i,:)=datas(i).ypoints2(firstInd-indBefore+1:firstInd+indAfter);
     
    ncon(i)=datas(i).numtrials1;
    nled(i)=datas(i).numtrials2;
end

% Make weighted average
scon=zeros(1,indBefore+indAfter);
sled=zeros(1,indBefore+indAfter);
for i=1:length(ncon)
    scon=scon+ncon(i)*allcon(i,:);
    sled=sled+nled(i)*allled(i,:);
end
scon=scon/sum(ncon);
sled=sled/sum(nled);

% Plot weighted average
figure(); 
plot(useX,scon,'Color','k');
hold on; 
plot(useX,sled,'Color','r');