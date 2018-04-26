function putTausTogetherFromAcrossFiles(xpoints,data1,data2,downSampFactor,alpha1,alpha2)
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

% baselineInds=find(xpoints>0 & xpoints<1);
% peakInds=find(xpoints>3.15 & xpoints<3.3);
% sharedBaseline=mean(mean(data1(:,baselineInds)));
% % Normalize by trial WITHOUT LED
% for i=1:size(data1,1)
%     height=mean(data1(i,peakInds))-mean(data1(i,baselineInds));
%     normFactor(i)=1/height;
% end
% % for i=1:size(data2,1)
% %     height=mean(data2(i,peakInds))-mean(data2(i,baselineInds));
% %     normFactor2(i)=1/height;
% % end
% for i=1:size(data1,1)
%     data1(i,:)=(data1(i,:)-mean(data1(i,baselineInds)))*normFactor(i);
%     data2(i,:)=(data2(i,:)-mean(data2(i,baselineInds)))*normFactor(i);
% end

baselineInds=find(xpoints>0 & xpoints<1);
peakInds=find(xpoints>3.15 & xpoints<3.3);
% Only for spontaneous activity analysis
% sharedBaseline=mean(mean(data1(:,baselineInds)));
% data1=data1+sharedBaseline;
% data2=data2+sharedBaseline;
% disp('sharedBaseline');
% disp(sharedBaseline);
% Normalize by trial WITHOUT LED
for i=1:size(data1,1)
    height=mean(data1(i,peakInds))-mean(data1(i,baselineInds));
    normFactor(i)=1/height;
%     normFactor(i)=1; % for spontaneous silencing
end
% for i=1:size(data2,1)
%     height=mean(data2(i,peakInds))-mean(data2(i,baselineInds));
%     normFactor2(i)=1/height;
% end
for i=1:size(data1,1)
    data1(i,:)=(data1(i,:)-mean(data1(i,baselineInds)))*normFactor(i);
    data2(i,:)=(data2(i,:)-mean(data2(i,baselineInds)))*normFactor(i);
end

% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(1,:)=data1(1,(xpoints>t1)&(xpoints<t2));
% dataY2(1,:)=data2(1,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(2,:)=data1(2,(xpoints>t1)&(xpoints<t2));
% dataY2(2,:)=data2(2,(xpoints>t1)&(xpoints<t2));
% t1=3.6-0.3;
% t2=t1+0.8;
% dataY1(3,:)=data1(3,(xpoints>t1)&(xpoints<t2));
% dataY2(3,:)=data2(3,(xpoints>t1)&(xpoints<t2));
% t1=3.9-0.3;
% t2=t1+0.8;
% dataY1(4,:)=data1(4,(xpoints>t1)&(xpoints<t2));
% dataY2(4,:)=data2(4,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(5,:)=data1(5,(xpoints>t1)&(xpoints<t2));
% dataY2(5,:)=data2(5,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(6,:)=data1(6,(xpoints>t1)&(xpoints<t2));
% dataY2(6,:)=data2(6,(xpoints>t1)&(xpoints<t2));
% t1=3.6-0.3;
% t2=t1+0.8;
% dataY1(7,:)=data1(7,(xpoints>t1)&(xpoints<t2));
% dataY2(7,:)=data2(7,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(8,:)=data1(8,(xpoints>t1)&(xpoints<t2));
% dataY2(8,:)=data2(8,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(9,:)=data1(9,(xpoints>t1)&(xpoints<t2));
% dataY2(9,:)=data2(9,(xpoints>t1)&(xpoints<t2));
% t1=3.6-0.3;
% t2=t1+0.8;
% dataY1(10,:)=data1(10,(xpoints>t1)&(xpoints<t2));
% dataY2(10,:)=data2(10,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(11,:)=data1(11,(xpoints>t1)&(xpoints<t2));
% dataY2(11,:)=data2(11,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(12,:)=data1(12,(xpoints>t1)&(xpoints<t2));
% dataY2(12,:)=data2(12,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(13,:)=data1(13,(xpoints>t1)&(xpoints<t2));
% dataY2(13,:)=data2(13,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(14,:)=data1(14,(xpoints>t1)&(xpoints<t2));
% dataY2(14,:)=data2(14,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(15,:)=data1(15,(xpoints>t1)&(xpoints<t2));
% dataY2(15,:)=data2(15,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;
% dataY1(16,:)=data1(16,(xpoints>t1)&(xpoints<t2));
% dataY2(16,:)=data2(16,(xpoints>t1)&(xpoints<t2));
% t1=3.6-0.3;
% t2=t1+0.8;
% dataY1(17,:)=data1(17,(xpoints>t1)&(xpoints<t2));
% dataY2(17,:)=data2(17,(xpoints>t1)&(xpoints<t2));
% t1=3.9-0.3;
% t2=t1+0.8;
% dataY1(18,:)=data1(18,(xpoints>t1)&(xpoints<t2));
% dataY2(18,:)=data2(18,(xpoints>t1)&(xpoints<t2));
% t1=3.3-0.3;
% t2=t1+0.8;

t1=3.2-0.3;
t2=t1+0.8;

dataY1(:,:)=data1(:,(xpoints>t1)&(xpoints<t2));
dataY2(:,:)=data2(:,(xpoints>t1)&(xpoints<t2));

xpoints=xpoints((xpoints>t1)&(xpoints<t2));
% trialDuration=0.8;
% 
% % if isempty(xpoints)
%     timestep=trialDuration/(size(dataY1,2));
%     xpoints=timestep/2:timestep:trialDuration;
% % end

% Use a larger bin size for PSTHs
for i=1:size(dataY1,1)
%     [dataY1{i},dataY2{i},outXpoints]=useLargerBinsize(dataY1{i},dataY2{i},xpoints,30);
    [newDataY1(i,:),newDataY2(i,:),newXpoints]=useLargerBinsize(dataY1(i,:),dataY2(i,:),xpoints,downSampFactor);
end

% Normalize PSTHs
% baselineInds=newXpoints>3.4 & newXpoints<3.5;
% for i=1:size(newDataY1,1)
%     newDataY1(i,:)=newDataY1(i,:)-mean(newDataY1(i,baselineInds));
%     newDataY2(i,:)=newDataY2(i,:)-mean(newDataY2(i,baselineInds));
% end

% baselineInds=newXpoints>3.4 & newXpoints<3.8;
% peakInds=newXpoints>3.2 & newXpoints<3.3;
% % Normalize by trial WITH LED
% for i=1:size(newDataY2,1)
%     height=mean(newDataY2(i,peakInds))-mean(newDataY2(i,baselineInds));
%     normFactor(i)=1/height;
% end
% for i=1:size(newDataY1,1)
%     newDataY1(i,:)=(newDataY1(i,:))*normFactor(i);
%     newDataY2(i,:)=(newDataY2(i,:))*normFactor(i);
% end


% Get S.E.M. measurements
stds1=std(newDataY1,1);
stds2=std(newDataY2,1);

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
errorbar(newXpoints,mean(newDataY1,1),stds1,'Color',[0.5 0.5 0.5]);
hold on;
errorbar(newXpoints,mean(newDataY2,1),stds2,'Color',[1 0.6 0.6]);
plot(newXpoints,mean(newDataY1,1),'Color','k');
plot(newXpoints,mean(newDataY2,1),'Color','r');
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