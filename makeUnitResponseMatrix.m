function makeUnitResponseMatrix(spikes,newassignsinfo)

norm=1;

a=unique(spikes.assigns);
s=unique(spikes.stimcond);

r=9;

[t1,t2,t3]=scriptForComparingMUA(filtspikes(spikes,0,'assigns',a(1),'stimcond',s(r)),[],[],[],0,5);
y1=zeros(length(a),size(t2,2));
y2=zeros(length(a),size(t2,2));
for i=1:length(a)
    [x,y1(i,:),y2(i,:)]=scriptForComparingMUA(filtspikes(spikes,0,'assigns',a(i),'stimcond',s(r)),[],[],[],0,5);
end

idx1=kmeans([y1 y2],10);

[s1,i1]=sort(idx1);
y1=y1(i1,:);
y2=y2(i1,:);

if norm==1
    m1=max(y1,[],2);
    m2=max(y2,[],2);
    y1(m1~=0,:)=y1(m1~=0,:)./repmat(m1(m1~=0),1,size(y1,2));
    y2(m2~=0,:)=y2(m2~=0,:)./repmat(m2(m2~=0),1,size(y2,2));
end

figure();
imagesc(y1);
title('LED Off');

figure();
imagesc(y2);
title('LED On');

% phasex=[1.2 1.388 1.563 1.738 1.913 2.138];
phasex=[4.038 4.237 4.413 4.588 4.737];

y2=y1;
allphaseslices=zeros(size(y2,1),1000);
allsets=zeros(size(y2,1)*length(phasex),1000);
for i=1:length(phasex)-1
    allphaseslices(:,1:sum(x>=phasex(i) & x<=phasex(i+1)+0.05))=allphaseslices(:,1:sum(x>=phasex(i) & x<=phasex(i+1)+0.05))+y2(:,x>=phasex(i) & x<=phasex(i+1)+0.05);
    allsets(i:length(phasex):end,1:sum(x>=phasex(i) & x<=phasex(i+1)+0.05))=y2(:,x>=phasex(i) & x<=phasex(i+1)+0.05);
end
figure();
imagesc(allphaseslices(:,1:20));
title('Period Spikes');

idx3=kmeans(allphaseslices,10);
[s3,i3]=sort(idx3);
allphaseslices=allphaseslices(i3,:);
for i=1:length(phasex)-1
    sub=allsets(i:length(phasex):end,:);
    allsets(i:length(phasex):end,:)=sub(i3,:);
end
figure();
imagesc(allphaseslices(:,1:20));
title('Period Spikes');
disp(newassignsinfo.calibrated_evCh(i3));

figure();
altogethersum=sum(allphaseslices(:,1:20),1);
plot(altogethersum);

figure(); 
imagesc(allsets(:,1:20));
title('Multiple Trials Phase');