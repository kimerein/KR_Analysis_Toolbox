function makeCSDtypeFig(dataMatrix,xtimes,oneWindow,zeroWindow)

normFirst=1;
notZeros=zeros(size(dataMatrix,1),1);
if normFirst==1
    for i=1:size(dataMatrix,1)
        currZero=mean(dataMatrix(i,xtimes>=zeroWindow(1) & xtimes<=zeroWindow(2)));
        dataMatrix(i,:)=dataMatrix(i,:)-currZero;
        currOne=mean(dataMatrix(i,xtimes>=oneWindow(1) & xtimes<=oneWindow(2)));
        if currOne~=0
            notZeros(i)=1;
        end
        dataMatrix(i,:)=dataMatrix(i,:)./currOne;
    end
end
subDataMatrix=dataMatrix(logical(notZeros),xtimes>=oneWindow(1) & xtimes<=zeroWindow(2));
dataMatrix=subDataMatrix(5:end,:);

h=(1/9)*ones(3);
Zsmooth=filter2(h,dataMatrix);
figure(); 
imagesc(Zsmooth);

% 2D linear interpolation
[x y]=size(Zsmooth);
y=1:y;
x=1:x;
[xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
dataInt=interp2(x,y,Zsmooth',xi,yi);

figure(); 
dataInt=dataInt';
imagesc(dataInt);
drawnow;
end

% csd=zeros(size(dataMatrix,1),size(dataMatrix,2));
% csd(1,:)=(2*dataMatrix(1,:)+dataMatrix(2,:))/3;
% csd(end,:)=(2*dataMatrix(end,:)+dataMatrix(end-1,:))/3;
% for i=2:size(dataMatrix,1)-1
%     csd(i,:)=(dataMatrix(i-1,:)+2*dataMatrix(i,:)+dataMatrix(i+1,:))/4;
% end
% for i=1:size(dataMatrix,1)
%     csd(i,:)=smooth(csd(i,:),binsize);
% end
% data1=csd(1:2:end,:);
% data2=csd(2:2:end,:);
% data1=diff(data1,2,1);
% data2=diff(data2,2,1);
% data=zeros(length(data1),2*size(data1,1));
% DataCounter=1;
% for i=1:size(data1,1)
%     data(:,DataCounter)=data1(i,:);
%     DataCounter=DataCounter+1;
%     if DataCounter>size(data2,1)
%         break
%     end
%     data(:,DataCounter)=data2(i,:);
%     DataCounter=DataCounter+1;
% end
% data=-data';
% 
% % 2D linear interpolation
% % [x y]=size(data);
% % y=1:y;
% % x=1:x;
% % [xi yi]=meshgrid(1:1:max(x),1:0.05:max(y));
% % dataInt=interp2(x,y,data',xi,yi);
% [x y]=size(data);
% y=1:y;
% x=1:x;
% [xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
% dataInt=interp2(x,y,data',xi,yi);
% 
% figure(); 
% dataInt=dataInt';
% imagesc(dataInt);
% drawnow;
% end