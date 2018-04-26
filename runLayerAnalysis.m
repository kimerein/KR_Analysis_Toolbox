function [subData2,sData2,delays,subData]=runLayerAnalysis()

% top=[3.28 3.3];
% bottom=[3.4 3.5];
% top=[1.18 1.2];
% bottom=[1.3 1.325];
top=[1.03 1.05];
bottom=[1.15 1.175];
% top=[1.32 1.34];
% bottom=[1.44 1.54];
% top=[1.28 1.3];
% bottom=[1.4 1.5];
cleanLayers=6:10;
smoothWindow=10;
nTimes=1;
smoothAcrossLs=0;
acrossLsSmoothWindow=2;

[dataY1,dataY2,xpoints]=loadInLayerData();
subData=normalizeFR_turningOff_acrossChs(dataY2,top,bottom);
sData=normalize_matrix_rows(xpoints,dataY2,top,bottom);

% ma=max(max(subData(cleanLayers,:)));
% mi=min(min(subData(cleanLayers,:)));
subData2=subData;
% subData2(subData2>ma)=ma;
% subData2(subData2<mi)=mi;
for j=1:nTimes
    for i=1:size(subData2,1)
        subData2(i,:)=smooth(subData2(i,:),smoothWindow);
    end
end
if smoothAcrossLs==1
    for i=1:size(subData,2)
        subData2(:,i)=smooth(subData2(:,i),acrossLsSmoothWindow);
    end
end
% [x y]=size(subData2);
% y=1:y;
% x=1:x;
% [xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
% interp_subData2=interp2(x,y,subData2',xi,yi);
figure(); 
temp=subData2;
ma=max(max(temp(cleanLayers,:)));
mi=min(min(temp(cleanLayers,:)));
temp(temp>ma)=ma;
temp(temp<mi)=mi;
imagesc(temp);
delays=getTimeToHalfMax(subData2,xpoints(2)-xpoints(1),0.5,top(2)-top(1)+0.003);

% ma=max(max(sData(cleanLayers,:)));
% mi=min(min(sData(cleanLayers,:)));
sData2=sData;
% sData2(sData2>ma)=ma;
% sData2(sData2<mi)=mi;
for j=1:nTimes
    for i=1:size(sData2,1)
%         sData2(i,:)=smooth(sData2(i,:),smoothWindow);
        newsData2(i,:)=downSampAv(sData2(i,:),3);
        newsData2(i,:)=smooth(newsData2(i,:),3);
    end
end
if smoothAcrossLs==1
    for i=1:size(sData2,2)
        sData2(:,i)=smooth(sData2(:,i),acrossLsSmoothWindow);
    end
end
[x y]=size(newsData2');
y=1:y;
x=1:x;
[xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
interp_sData2=interp2(x,y,newsData2,xi,yi);
for i=1:size(interp_sData2,1)
    interp_sData2(i,:)=smooth(interp_sData2(i,:),3);
end
% for i=1:size(interp_sData2,2)
%     interp_sData2(:,i)=smooth(interp_sData2(:,i),3);
% end
figure(); 
imagesc(interp_sData2);
% imagesc(sData2);
figure(); 
temp=sData2;
ma=max(max(temp(cleanLayers,:)));
mi=min(min(temp(cleanLayers,:)));
temp(temp>ma)=ma;
temp(temp<mi)=mi;
imagesc(temp);

figure();
plot(1:16,delays);
