function subData=normalizeFR_turningOff_acrossChs(data,top,bottom)

% startTop=1.1;
% endTop=1.2;
% startTop=1.1;
% endTop=1.3;

startTop=top(1);
endTop=top(2);

% startBottom=1.28;
% endBottom=1.31;
% startBottom=1.38;
% endBottom=1.41;
% startBottom=1.7;
% endBottom=1.8;

startBottom=bottom(1);
endBottom=bottom(2);

smoothWindow=1;
smoothAcrossLs=0;
smoothLsWindow=2;

startTop=startTop/(5/size(data,2));
endTop=endTop/(5/size(data,2));
startBottom=startBottom/(5/size(data,2));
endBottom=endBottom/(5/size(data,2));
startTop=floor(startTop);
endTop=floor(endTop);
startBottom=floor(startBottom);
endBottom=floor(endBottom);

for i=1:16
    data(i,:)=smooth(data(i,:),smoothWindow);
end
if smoothAcrossLs==1
    for i=1:size(data,2)
        data(:,i)=smooth(data(:,i),2);
    end
end

subData=data(:,startTop:endBottom);

startTop2=startTop-startTop+1;
endTop2=endTop-startTop+1;
startBottom2=startBottom-startTop+1;
endBottom2=endBottom-startTop+1;
startTop=startTop2;
endTop=endTop2;
startBottom=startBottom2;
endBottom=endBottom2;

for i=1:16
    postAvs(i)=mean(subData(i,startBottom:endBottom));
end

for i=1:16
    subData(i,:)=subData(i,:)-postAvs(i);
end
subData(subData<0)=0;

for i=1:16
    preAvs(i)=mean(subData(i,startTop:endTop));
end

for i=1:16
    subData(i,:)=(200/preAvs(i))*subData(i,:);
end

for i=1:16
    preAvs(i)=mean(subData(i,startTop:endTop));
    postAvs(i)=mean(subData(i,startBottom:endBottom));
end
% disp(preAvs);
% disp(postAvs);

% figure();
% imagesc(subData);
