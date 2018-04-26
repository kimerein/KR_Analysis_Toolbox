function [subData2,sData2,delays,subData]=runLayerAnalysis2(dataY2,xpoints)

% top=[3.28 3.3];
% bottom=[3.4 3.5];
% top=[1.18 1.2];
% bottom=[1.3 1.325];
top=[1.28 1.298];
bottom=[1.5 1.6];
% top=[1.32 1.34];
% bottom=[1.44 1.54];
% top=[1.28 1.3];
% bottom=[1.4 1.5];
nTimes=1;

if isempty(xpoints)
    [dataY1,dataY2,xpoints]=loadInLayerData();
end
subData=normalizeFR_turningOff_acrossChs(dataY2,top,bottom);
sData=normalize_matrix_rows(xpoints,dataY2,top,bottom);

subData2=subData;
timeW=size(subData2,2)*0.001;
for j=1:nTimes
    for i=1:size(subData2,1)
        newsubData2(i,:)=downSampAv(subData2(i,:),6);
    end
end
[x y]=size(newsubData2');
y=1:y;
x=1:x;
[xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
interp_subData2=interp2(x,y,newsubData2,xi,yi);
figure(); 
imagesc(interp_subData2);
delays=getTimeToHalfMax(interp_subData2,timeW/size(interp_subData2,2),0.5,top(2)-top(1)+0.003);

sData2=sData;
for j=1:nTimes
    for i=1:size(sData2,1)
        newsData2(i,:)=downSampAv(sData2(i,:),6);
    end
end
[x y]=size(newsData2');
y=1:y;
x=1:x;
[xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
interp_sData2=interp2(x,y,newsData2,xi,yi);
figure(); 
imagesc(interp_sData2(50:end,:));

figure();
plot(1:16,delays);
