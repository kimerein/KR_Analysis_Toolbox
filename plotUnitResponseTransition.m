function [togethery,x1]=plotUnitResponseTransition(spikes,ass,binsize,grp2_stimconds)

% fileInd_grp1=[71:83];
% fileInd_grp1=[59:69];
% fileInd_grp1=[29:51];
fileInd_grp1=[69:89];
fileInd_grp2=[52:141];
% fileInd_grp2=[70:148];
% fileInd_grp2=[120:176];

spikes1=filtspikes(spikes,0,'fileInd',fileInd_grp1);
ys=nan(length(1:12),length(1:10000));
for i=1:12
    [~,~,~,x1_temp,y1_temp,~]=psth_wStdev_valuesOnly(filtspikes(spikes1,0,'assigns',ass,'stimcond',i),binsize);
    ys(i,1:length(y1_temp))=y1_temp;
end
firstNan=find(isnan(ys(1,:)),1,'first');
ys=ys(:,1:firstNan-1);
[~,~,~,x1_1,y1_1,~]=psth_wStdev_valuesOnly(filtspikes(spikes1,0,'assigns',ass,'stimcond',[1:4]),binsize);
[~,~,~,x2_1,y2_1,~]=psth_wStdev_valuesOnly(filtspikes(spikes1,0,'assigns',ass,'stimcond',[5:8]),binsize);
[~,~,~,x3_1,y3_1,~]=psth_wStdev_valuesOnly(filtspikes(spikes1,0,'assigns',ass,'stimcond',[9:12]),binsize);
% [~,~,~,x1_1,y1_1,~]=psth_wStdev_valuesOnly(filtspikes(spikes1,0,'assigns',ass,'stimcond',[1 4 7 10]),binsize);
% [~,~,~,x2_1,y2_1,~]=psth_wStdev_valuesOnly(filtspikes(spikes1,0,'assigns',ass,'stimcond',[2 5 8 11]),binsize);
% [~,~,~,x3_1,y3_1,~]=psth_wStdev_valuesOnly(filtspikes(spikes1,0,'assigns',ass,'stimcond',[3 6 9 12]),binsize);

spikes2=filtspikes(spikes,0,'fileInd',fileInd_grp2);
% [~,~,~,x1,y1,~]=psth_wStdev_valuesOnly(filtspikes(spikes2,0,'assigns',ass,'stimcond',[1 4 7 10]),binsize);
% [~,~,~,x2,y2,~]=psth_wStdev_valuesOnly(filtspikes(spikes2,0,'assigns',ass,'stimcond',[2 5 8 11]),binsize);
% [~,~,~,x3,y3,~]=psth_wStdev_valuesOnly(filtspikes(spikes2,0,'assigns',ass,'stimcond',[3 6 9 12]),binsize);
[~,~,~,x1,y1,~]=psth_wStdev_valuesOnly(filtspikes(spikes2,0,'assigns',ass,'stimcond',grp2_stimconds(1)),binsize);
[~,~,~,x2,y2,~]=psth_wStdev_valuesOnly(filtspikes(spikes2,0,'assigns',ass,'stimcond',grp2_stimconds(2)),binsize);
[~,~,~,x3,y3,~]=psth_wStdev_valuesOnly(filtspikes(spikes2,0,'assigns',ass,'stimcond',grp2_stimconds(3)),binsize);

fig_min=min([y1_1 y2_1 y3_1 y1 y2 y3]);
fig_max=max([y1_1 y2_1 y3_1 y1 y2 y3]);

figure();
plot(x1_1(x1_1>=0.3 & x1_1<=1.8)-0.3,y1_1(x1_1>=0.3 & x1_1<=1.8),'Color','b');
hold on;
plot(x2_1(x2_1>=0.3 & x2_1<=1.8)-0.3+1.5,y2_1(x2_1>=0.3 & x2_1<=1.8),'Color','b');
plot(x3_1(x3_1>=0.3 & x3_1<=1.8)-0.3+3,y3_1(x3_1>=0.3 & x3_1<=1.8),'Color','b');
% ylim([fig_min fig_max]);

m{1}=mean(ys(1:4,x1_temp>=0.3 & x1_temp<=1.8),2);
m{2}=mean(ys(5:8,x1_temp>=0.3 & x1_temp<=1.8),2);
m{3}=mean(ys(9:12,x1_temp>=0.3 & x1_temp<=1.8),2);
[~,mi]=max([mean(m{1}) mean(m{2}) mean(m{3})]);
disp(mi);
i=[1 2 3];
other=i(~ismember(i,mi));
[~,p1]=ttest(m{mi},m{other(1)});
[~,p2]=ttest(m{mi},m{other(2)});
figure();
for i=1:length(m)
    if any([p1 p2]<0.05) && i==mi
        scatter(ones(size(m{i})).*i,m{i},[],'r');
    else
        scatter(ones(size(m{i})).*i,m{i},[],'b');
    end
    hold on;
end
xlim([0.5 3.5]);
return
figure(); 
plot(x1,y1,'Color','k');
hold on; 
plot(x2(x2>=0 & x2<=3)+1.5,y2(x2>=0 & x2<=3),'Color','r');
plot(x2(x2>=3)-3,y2(x2>=3),'Color','r');
plot(x3(x3>=0 & x3<=1.5)+3,y3(x3>=0 & x3<=1.5),'Color','g');
plot(x3(x3>=1.5 & x3<=3)-1.5,y3(x3>=1.5 & x3<=3),'Color','g');
plot(x3(x3>=3)-1.5,y3(x3>=3),'Color','g');
togethery=nan(3,length(x1(x1>=0 & x1<=1.5))*3);
x1=x1(1):x1(2)-x1(1):4.5;
togethery(1,x1>=0 & x1<=3.5)=y1;
togethery(2,x1>=0 & x1<=0.5)=y2(x2>=3);
togethery(2,x1>=1.5 & x1<=4.5)=y2(x2>=0 & x2<=3);
togethery(3,x1>=0 & x1<=1.5)=y3(x3>=1.5 & x3<=3);
togethery(3,x1>=1.5 & x1<=2)=y3(x3>=3);
togethery(3,x1>=3 & x1<=4.5)=y3(x3>=0 & x3<=1.5);
% ylim([fig_min fig_max]);