function [newXpoints,newDataY1,stds1,newDataY2,stds2]=putPSTHsTogetherFromAcrossFiles(xpoints,dataY1,dataY2,downSampFactor,alpha1,alpha2)
% Use these data sets
% xpoints=xpoints1;

% dataY1{1}=ypoints1_120211;
% dataY1{2}=ypoints1_120221;
% dataY1{3}=ypoints1_120228;
% dataY1{4}=ypoints1_120319;
% dataY1{5}=ypoints1_120321;
% 
% dataY2{1}=ypoints2_120211;
% dataY2{2}=ypoints2_120221;
% dataY2{3}=ypoints2_120228;
% dataY2{4}=ypoints2_120319;
% dataY2{5}=ypoints2_120321;

% dataY1(1,:)=ypoints1_120211;
% dataY1(2,:)=ypoints1_120221;
% dataY1(3,:)=ypoints1_120228;
% dataY1(4,:)=ypoints1_120319;
% dataY1(5,:)=ypoints1_120321;
% 
% dataY2(1,:)=ypoints2_120211;
% dataY2(2,:)=ypoints2_120221;
% dataY2(3,:)=ypoints2_120228;
% dataY2(4,:)=ypoints2_120319;
% dataY2(5,:)=ypoints2_120321;

trialDuration=5;

if isempty(xpoints)
    timestep=trialDuration/(size(dataY1,2));
    xpoints=timestep/2:timestep:trialDuration;
end

% Use a larger bin size for PSTHs
for i=1:size(dataY1,1)
%     [dataY1{i},dataY2{i},outXpoints]=useLargerBinsize(dataY1{i},dataY2{i},xpoints,30);
    [newDataY1(i,:),newDataY2(i,:),newXpoints]=useLargerBinsize(dataY1(i,:),dataY2(i,:),xpoints,downSampFactor);
end

% Normalize PSTHs
baselineInds=find(newXpoints>0 & newXpoints<0.48);
% baselineInds=find(newXpoints>0.85 & newXpoints<1);
% baselineInds=find(newXpoints>0 & newXpoints<0.4);

% baselineInds=find(newXpoints>2 & newXpoints<3);
% peakInds=find(newXpoints>1.15 & newXpoints<1.78);
% peakInds=find(newXpoints>1.05 & newXpoints<1.28);
peakInds=find(newXpoints>1.2 & newXpoints<1.3);
% peakInds{2}=find(newXpoints>3.07 & newXpoints<3.13);
% peakInds{1}=find(newXpoints>3.2 & newXpoints<3.6);
% peakInds=find(newXpoints>3 & newXpoints<3.1);
% peakInds=find(newXpoints>1.05 & newXpoints<1.15);
sharedBaseline=mean(mean(newDataY1(:,baselineInds)));
% sharedBaseline=mean(mean(newDataY1(1,baselineInds)));
% Normalize by trial WITHOUT LED
for i=1:size(newDataY1,1)
%     height=mean(newDataY1(i,peakInds))-mean(newDataY1(i,baselineInds));
    height=mean(newDataY2(i,peakInds))-mean(newDataY2(i,baselineInds));
    normFactor(i)=1/height;
%     normFactor(i)=1; % for spontaneous activity study
end
% for i=1:size(newDataY2,1)
%     height=mean(newDataY2(i,peakInds))-mean(newDataY2(i,baselineInds));
%     normFactor2(i)=1/height;
% %     normFactor(i)=1; % for spontaneous activity study
% end
for i=1:size(newDataY1,1)
    newDataY1(i,:)=(newDataY1(i,:)-mean(newDataY1(i,baselineInds)))*normFactor(i);
    newDataY2(i,:)=(newDataY2(i,:)-mean(newDataY2(i,baselineInds)))*normFactor(i);
end
newDataY1=newDataY1+sharedBaseline;
newDataY2=newDataY2+sharedBaseline;

% Get S.E.M. measurements
stds1=std(newDataY1,[],1);
stds2=std(newDataY2,[],1);

% All lines figures
figure();
cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
    [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
    [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2]};
for i=1:size(newDataY1,1)
    plot(newXpoints,newDataY1(i,:),'Color',cs{i});
    hold on;
end
figure();
for i=1:size(newDataY2,1)
    plot(newXpoints,newDataY2(i,:),'Color',cs{i});
    hold on;
end

% S.E.M. figure
figure();
% errorbar(newXpoints,mean(newDataY1,1),stds1,'Color',[0.5 0.5 0.5]);
hax=axes();
hl=plot(newXpoints,mean(newDataY1,1),'Color','k');
addErrBar(newXpoints,mean(newDataY1,1),stds1,'y',hax,hl);
% errorbar(newXpoints,mean(newDataY1,1),stds1,'Color',[0.5 0.5 0.5]);
hold on;
hl=plot(newXpoints,mean(newDataY2,1),'Color','r');
addErrBar(newXpoints,mean(newDataY2,1),stds2,'y',hax,hl);
% errorbar(newXpoints,mean(newDataY2,1),stds2,'Color',[1 0.6 0.6]);
% plot(newXpoints,mean(newDataY1,1),'Color','k');
% plot(newXpoints,mean(newDataY2,1),'Color','r');
%hold off;

pvals=zeros(length(newXpoints));
for i=1:length(newXpoints)
    [h,p,ci]=ttest(newDataY1(:,i),newDataY2(:,i));
    pvals(i)=p;
end

maxBoth=max(max(mean(newDataY1,1)),max(mean(newDataY2,1)));
minBoth=min(min(mean(newDataY1,1)),min(mean(newDataY2,1)));
halfX=(newXpoints(2)-newXpoints(1))/2;
for i=1:length(newXpoints)
    if pvals(i)<alpha1
        line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth maxBoth],'Color',[0.5 0.5 0.5],'LineWidth',3);
    end
    above=(maxBoth-minBoth)/15;
    if pvals(i)<alpha2
        line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth+above maxBoth+above],'Color','k','LineWidth',2);
    end
end