function [returnX,returnLED,returnCon]=averageTausAcrossExpts(datas)

useSubtraction=0;
weighted=1;

indBefore=500;
indAfter=120;
% indBefore=1000;
% indAfter=1000;

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
    if useSubtraction==1
        allled(i,:)=mean(allcon(i,:)).*ones(size(allcon(i,:)))-(allcon(i,:)-allled(i,:));
    end      
    ncon(i)=datas(i).numtrials1;
    nled(i)=datas(i).numtrials2;
end

% Make weighted average
scon=zeros(1,indBefore+indAfter);
sled=zeros(1,indBefore+indAfter);
for i=1:length(ncon)
    if weighted==1
        scon=scon+ncon(i)*allcon(i,:);
        sled=sled+nled(i)*allled(i,:);
    else
        scon=scon+allcon(i,:);
        sled=sled+allled(i,:);
    end
end
scon=scon/sum(ncon);
sled=sled/sum(nled);
returnLED=sled;
returnX=useX;
returnCon=scon;

% Plot weighted average
figure(); 
plot(useX,scon,'Color','k');
hold on; 
plot(useX,sled,'Color','r');

% Down Sample
for i=1:length(ncon)
    newcon(i,:)=downSampAv(allcon(i,:),4);
    newled(i,:)=downSampAv(allled(i,:),4);
end
newX=downSampAv(useX,4);
clear allcon allled useX
useX=newX;
allcon=newcon;
allled=newled;

baseAv=0;
for i=1:length(ncon)
    currbase=mean(allcon(i,ismember(useX,-500:1:-100)));
%     allcon(i,:)=allcon(i,:)-currbase;
    allcon(i,:)=allcon(i,:)*(60/currbase);
    baseAv=baseAv+currbase;
%     allled(i,:)=allled(i,:)-currbase;
    allled(i,:)=allled(i,:)*(60/currbase);
end
% baseAv=baseAv/length(ncon);
% for i=1:length(ncon)
%     allcon(i,:)=allcon(i,:)+baseAv;
%     allled(i,:)=allled(i,:)+baseAv;
% end


% Make unweighted average with error
scon=zeros(1,size(allcon,2));
sled=zeros(1,size(allcon,2));
for i=1:length(ncon)
    scon=scon+allcon(i,:);
    sled=sled+allled(i,:);
end
scon=scon./size(allcon,1);
sled=sled./size(allled,1);
scon=mean(allcon,1);
sled=mean(allled,1);
scon_std=std(allcon,[],1);
sled_std=std(allled,[],1);
    
% Plot weighted average
figure(); 
hax=axes();
hl=plot(useX,scon,'Color','k');
addErrBar(useX,scon,scon_std,'y',hax,hl);
hold on; 
hl=plot(useX,sled,'Color','r');
addErrBar(useX,sled,sled_std,'y',hax,hl);