function [dataY1,dataY2,xpoints]=loadInLayerData()

dataDir='W:\Analysis Computer\Layers MUA\Awake\Mawake49\';
fname='led0vs5';

for i=1:4
    for j=1:4
        fullname=[dataDir 'T' num2str(i) '_evch' num2str(j) '_' fname '.mat'];
        a=load(fullname);
        dataY1((i-1)*4+j,:)=a.data.ypoints1;
        dataY2((i-1)*4+j,:)=a.data.ypoints2;
    end
end
xpoints=a.data.xpoints;
