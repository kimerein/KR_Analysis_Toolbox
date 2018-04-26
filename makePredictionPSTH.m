function togethery=makePredictionPSTH(spikes,prefstimset,binsize,doGp1)

a=unique(spikes.assigns);

if doGp1==1
    [~,~,~,x,psths(1,1,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==1),'stimcond',1:4),binsize,0);
    [~,~,~,x,psths(1,2,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==1),'stimcond',5:8),binsize,0);
    [~,~,~,x,psths(1,3,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==1),'stimcond',9:12),binsize,0);
    [~,~,~,x,psths(2,1,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==2),'stimcond',1:4),binsize,0);
    [~,~,~,x,psths(2,2,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==2),'stimcond',5:8),binsize,0);
    [~,~,~,x,psths(2,3,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==2),'stimcond',9:12),binsize,0);
    [~,~,~,x,psths(3,1,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==3),'stimcond',1:4),binsize,0);
    [~,~,~,x,psths(3,2,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==3),'stimcond',5:8),binsize,0);
    [~,~,~,x,psths(3,3,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==3),'stimcond',9:12),binsize,0);
    
    figure();
    plot(linspace(0,1.8,size(psths,3)),reshape(psths(1,1,:),1,length(psths(1,1,:))),'Color','k');
    hold on;
    plot(linspace(0,1.8,size(psths,3)),reshape(psths(2,1,:),1,length(psths(1,1,:))),'Color','r');
    plot(linspace(0,1.8,size(psths,3)),reshape(psths(3,1,:),1,length(psths(1,1,:))),'Color','g');
    title('Stimcond 1');
    
    figure();
    plot(linspace(0,1.8,size(psths,3)),reshape(psths(1,2,:),1,length(psths(1,1,:))),'Color','k');
    hold on;
    plot(linspace(0,1.8,size(psths,3)),reshape(psths(2,2,:),1,length(psths(1,1,:))),'Color','r');
    plot(linspace(0,1.8,size(psths,3)),reshape(psths(3,2,:),1,length(psths(1,1,:))),'Color','g');
    title('Stimcond 2');
    
    figure();
    plot(linspace(0,1.8,size(psths,3)),reshape(psths(1,3,:),1,length(psths(1,1,:))),'Color','k');
    hold on;
    plot(linspace(0,1.8,size(psths,3)),reshape(psths(2,3,:),1,length(psths(1,1,:))),'Color','r');
    plot(linspace(0,1.8,size(psths,3)),reshape(psths(3,3,:),1,length(psths(1,1,:))),'Color','g');
    title('Stimcond 3');
    togethery=psths;
else
    [~,~,~,x,psths(1,1,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==1),'stimcond',[1 5 9 13]),binsize,0);
    [~,~,~,x,psths(1,2,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==1),'stimcond',1+[1 5 9 13]),binsize,0);
    [~,~,~,x,psths(1,3,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==1),'stimcond',2+[1 5 9 13]),binsize,0);
    [~,~,~,x,psths(2,1,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==2),'stimcond',[1 5 9 13]),binsize,0);
    [~,~,~,x,psths(2,2,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==2),'stimcond',1+[1 5 9 13]),binsize,0);
    [~,~,~,x,psths(2,3,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==2),'stimcond',2+[1 5 9 13]),binsize,0);
    [~,~,~,x,psths(3,1,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==3),'stimcond',[1 5 9 13]),binsize,0);
    [~,~,~,x,psths(3,2,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==3),'stimcond',1+[1 5 9 13]),binsize,0);
    [~,~,~,x,psths(3,3,:)]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'assigns',a(prefstimset==3),'stimcond',2+[1 5 9 13]),binsize,0);
    
    x1=x(1):x(2)-x(1):4.5;
    y=nan(3,length(x(x>=0 & x<=1.5))*3);
    togethery=nan(3,length(x(x>=0 & x<=1.5))*3);
    for i=1:3
        temp=makeVecsSameLength(y(1,x1>=0 & x1<=3.5),psths(i,1,:));
        y(1,x1>=0 & x1<=3.5)=temp;
        temp=makeVecsSameLength(y(2,x1>=0 & x1<=0.5),psths(i,2,x>=3));
        y(2,x1>=0 & x1<=0.5)=temp;
        temp=makeVecsSameLength(y(2,find(x1>=1.5,1,'first'):end),psths(i,2,x>=0 & x<=3));
        y(2,find(x1>=1.5,1,'first'):end)=temp;
        temp=makeVecsSameLength(y(3,x1>=0 & x1<=2),psths(i,3,x>=1.5 & x<=3.5));
        y(3,x1>=0 & x1<=2)=temp;
        temp=makeVecsSameLength(y(3,find(x1>=3,1,'first'):end),psths(i,3,x>=0 & x<=1.5));
        y(3,find(x1>=3,1,'first'):end)=temp;
        togethery(i,:)=nanmean(y,1);
    end
    
    figure();
    plot(x1,togethery(1,:),'Color','k');
    hold on;
    plot(x1,togethery(2,:),'Color','r');
    plot(x1,togethery(3,:),'Color','g');
end

end

function temp=makeVecsSameLength(a,b)

% Makes vec b same length as vec a

temp=nan(1,length(a));
if length(a)>length(b)
    temp(1:length(b))=b;
elseif length(b)>length(a)
    temp=b(1:length(a));
else
    temp=b;
end

end