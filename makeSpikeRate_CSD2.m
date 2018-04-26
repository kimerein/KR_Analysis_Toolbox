function dataInt=makeSpikeRate_CSD2(data)

csd=zeros(size(data,1),size(data,2));
csd(1,:)=(2*data(1,:)+data(2,:))/3;
csd(16,:)=(2*data(16,:)+data(15,:))/3;
for i=2:15
    csd(i,:)=(data(i-1,:)+2*data(i,:)+data(i+1,:))/4;
end
for i=1:16
%     csd(i,:)=smooth(csd(i,:),3);
    csd(i,:)=smooth(csd(i,:),5);
end
data1=csd(1:2:16,:);
data2=csd(2:2:16,:);

% for i=1:size(data1,1)-1
%     data1temp(i,:)=mean([data1(i,:); data1(i+1,:)],1);
% end
% data1=data1temp;
% for i=1:size(data1,1)-1
%     data1temp(i,:)=mean([data1(i,:); data1(i+1,:)],1);
% end
% data1=data1temp;

% for i=1:size(data2,1)-1
%     data2temp(i,:)=mean([data2(i,:); data2(i+1,:)],1);
% end
% data2=data2temp;
% for i=1:size(data2,1)-1
%     data2temp(i,:)=mean([data2(i,:); data2(i+1,:)],1);
% end
% data2=data2temp;

data=zeros(length(data1),2*size(data1,1));
DataCounter=1;
for i=1:size(data1,1)
    data(:,DataCounter)=data1(i,:);
    DataCounter=DataCounter+1;
    data(:,DataCounter)=data2(i,:);
    DataCounter=DataCounter+1;
end
% data=-data';
data=data';

[x y]=size(data);
y=1:y;
x=1:x;
[xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
dataInt=interp2(x,y,data',xi,yi);

figure();
dataInt=dataInt';
imagesc(dataInt);
drawnow;