function [x_new,t1_new,t2_new,t3_new]=plotPvaluePSTH(x,t1,t2,t3,downSampBin)

doWhen=0; % instead of doing rate
norm=1;
normToMean=1;

for i=1:size(t1,1)
    t1_new(i,:)=downSampAv(t1(i,:),downSampBin);
    t2_new(i,:)=downSampAv(t2(i,:),downSampBin);
    t3_new(i,:)=downSampAv(t3(i,:),downSampBin);
    if norm==1
        t1_new(i,:)=t1_new(i,:)-min(t1_new(i,:));
        t1_new(i,:)=t1_new(i,:)./max(t1_new(i,:));
        t2_new(i,:)=t2_new(i,:)-min(t2_new(i,:));
        t2_new(i,:)=t2_new(i,:)./max(t2_new(i,:));
        t3_new(i,:)=t3_new(i,:)-min(t3_new(i,:));
        t3_new(i,:)=t3_new(i,:)./max(t3_new(i,:));
    end
end
x_new=downSampAv(x,downSampBin);

if normToMean==1
    t1_new(:,x_new>=0 & x_new<=1.5)=t1_new(:,x_new>=0 & x_new<=1.5)-nanmean(nanmean(t1_new(:,x_new>=0.5 & x_new<=1.5),2),1);
    t2_new(:,x_new>=0 & x_new<=1.5)=t2_new(:,x_new>=0 & x_new<=1.5)-nanmean(nanmean(t2_new(:,x_new>=0.5 & x_new<=1.5),2),1);
    t3_new(:,x_new>=0 & x_new<=1.5)=t3_new(:,x_new>=0 & x_new<=1.5)-nanmean(nanmean(t3_new(:,x_new>=0.5 & x_new<=1.5),2),1);
    t1_new(:,x_new>=1.5 & x_new<=3)=t1_new(:,x_new>=1.5 & x_new<=3)-nanmean(nanmean(t1_new(:,x_new>=2 & x_new<=3),2),1);
    t2_new(:,x_new>=1.5 & x_new<=3)=t2_new(:,x_new>=1.5 & x_new<=3)-nanmean(nanmean(t2_new(:,x_new>=2 & x_new<=3),2),1);
    t3_new(:,x_new>=1.5 & x_new<=3)=t3_new(:,x_new>=1.5 & x_new<=3)-nanmean(nanmean(t3_new(:,x_new>=2 & x_new<=3),2),1);
    t1_new(:,x_new>=3 & x_new<=4.5)=t1_new(:,x_new>=3 & x_new<=4.5)-nanmean(nanmean(t1_new(:,x_new>=3.5 & x_new<=4.5),2),1);
    t2_new(:,x_new>=3 & x_new<=4.5)=t2_new(:,x_new>=3 & x_new<=4.5)-nanmean(nanmean(t2_new(:,x_new>=3.5 & x_new<=4.5),2),1);
    t3_new(:,x_new>=3 & x_new<=4.5)=t3_new(:,x_new>=3 & x_new<=4.5)-nanmean(nanmean(t3_new(:,x_new>=3.5 & x_new<=4.5),2),1);
end

figure(); 
plot(x_new,t1_new,'Color','k');
hold on; 
plot(x_new,t2_new,'Color','r');
plot(x_new,t3_new,'Color','g');


figure(); 
hax=axes();
hl=plot(x_new,nanmean(t1_new,1),'Color','k');
hold on; 
hl=plot(x_new,nanmean(t2_new,1),'Color','r');
hl=plot(x_new,nanmean(t3_new,1),'Color','g');

figure(); 
hax=axes();
hl=plot(x_new,nanmean(t1_new,1),'Color','k');
addErrBar(x_new,nanmean(t1_new,1),nanstd(t1_new,[],1)./sqrt(size(t1_new,1)),'y',hax,hl);
hold on; 
hl=plot(x_new,nanmean(t2_new,1),'Color','r');
addErrBar(x_new,nanmean(t2_new,1),nanstd(t2_new,[],1)./sqrt(size(t2_new,1)),'y',hax,hl);
hl=plot(x_new,nanmean(t3_new,1),'Color','g');
addErrBar(x_new,nanmean(t3_new,1),nanstd(t3_new,[],1)./sqrt(size(t3_new,1)),'y',hax,hl);


pvalPSTH=ones(3,size(t1_new,2));
v=ones(3,size(t1_new,2));
xbin=x_new(2)-x_new(1);
% disp(size(t1_new,2));
for i=1:size(t1_new,2)
%     disp(i);
    if doWhen==0
%         [~,p]=ttest(t1_new(:,i),t2_new(:,i));
        a=t1_new(:,i);
        b=t2_new(:,i);
        v(1,i)=std(a(~isnan(a)));
        p=ranksum(a(~isnan(a)),b(~isnan(b)));
        if ~isnan(p)
            pvalPSTH(1,i)=p;
        else
            disp('stop');
        end
%         [~,p]=ttest(t2_new(:,i),t3_new(:,i));
        a=t2_new(:,i);
        b=t3_new(:,i);
        v(2,i)=std(a(~isnan(a)));
        v(3,i)=std(b(~isnan(b)));
        p=ranksum(a(~isnan(a)),b(~isnan(b)));
        if ~isnan(p)
            pvalPSTH(2,i)=p;
        end
%         [~,p]=ttest(t1_new(:,i),t3_new(:,i));
        a=t1_new(:,i);
        b=t3_new(:,i);
        p=ranksum(a(~isnan(a)),b(~isnan(b)));
        if ~isnan(p)
            pvalPSTH(3,i)=p;
        end
    else
        xInds=find(x>=(x_new(i)-xbin/2) & x<=(x_new(i)+xbin/2));
        spiketimes1=[];
        for j=1:length(xInds)
            spiketimes1=[spiketimes1 ones(1,floor(t1(xInds(j))*1)).*x(xInds(j))];
        end
        spiketimes2=[];
        for j=1:length(xInds)
            spiketimes2=[spiketimes2 ones(1,floor(t2(xInds(j))*1)).*x(xInds(j))];
        end
        spiketimes3=[];
        for j=1:length(xInds)
            spiketimes3=[spiketimes3 ones(1,floor(t3(xInds(j))*1)).*x(xInds(j))];
        end
        if isempty(spiketimes1) || isempty(spiketimes2)
            p=1;
        else
            p=ranksum(spiketimes1,spiketimes2);
        end
        if p==0
            disp('stop');
        end
        if ~isnan(p)
            pvalPSTH(1,i)=p;
        end
        if isempty(spiketimes2) || isempty(spiketimes3)
            p=1;
        else
            p=ranksum(spiketimes2,spiketimes3);
        end
        if ~isnan(p)
            pvalPSTH(2,i)=p;
        end
        if isempty(spiketimes1) || isempty(spiketimes3)
            p=1;
        else
            p=ranksum(spiketimes1,spiketimes3);
        end
        if ~isnan(p)
            pvalPSTH(3,i)=p;
        end
    end
end

% figure();
% plot(downSampAv(x_new,1),downSampAv(v(1,:),1),'Color','k');
% figure();
% plot(downSampAv(x_new,1),downSampAv(v(2,:),1),'Color','r');
% figure();
% plot(downSampAv(x_new,1),downSampAv(v(3,:),1),'Color','g');

figure();
plot(downSampAv(x_new,1),downSampAv(pvalPSTH(1,:),1),'Color','k');
figure();
plot(downSampAv(x_new,1),downSampAv(pvalPSTH(2,:),1),'Color','r');
figure();
plot(downSampAv(x_new,1),downSampAv(pvalPSTH(3,:),1),'Color','g');