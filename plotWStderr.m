function [data1,data2,isBigEnough]=plotWStderr(data1,data2,trialDuration,c1,c2)


norm=1;
normWindow=[1 3];
% normWindow=[3.5 4];
% normWindow=[4 6.5];
baseSubtract=1;
% baseWindow=[0 2.5];
% baseWindow=[0 3];
% baseWindow=[3.9 4];
baseWindow=[0 0.7];
% baseWindow=[3.7 3.9];
% baseWindow=[3 3.7];
% baseWindow=[0 3];
% baseWindow=[2.25 3.5];
% ds=[];
% ds=2;
% ds=3;
% ds=20;
ds=1;
doMed=0;
doFill=1;

if length(trialDuration)>1
    t=trialDuration;
else
%     t=linspace(0,trialDuration,size(data1,2));
%     t=linspace(1,trialDuration-1,size(data1,2));
    t=linspace(0.5,trialDuration-0.5,size(data1,2));
%     t=linspace(0.0825,trialDuration-0.0825,size(data1,2));
end


% % Test whether response is sustained
% useUnit=zeros(1,size(data2,1));
% for i=1:size(data2,1)
%     currBase=data2(i,t>=0 & t<=3);
%     currResponse=data2(i,t>=5 & t<=6);
%     p=ranksum(currBase,currResponse);
%     if p<0.05 && nanmean(currBase)*1<nanmean(currResponse)
%         useUnit(i)=1;
%     end
% end
% data1=data1(useUnit==1,:);
% data2=data2(useUnit==1,:);


if baseSubtract==1
    for i=1:size(data1,1)
        base1=nanmean(data1(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
        base2=nanmean(data2(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
        data1(i,:)=data1(i,:)-base1;
        data2(i,:)=data2(i,:)-base2;
    end
end

% % thresh=0.15*10^4;
% % thresh=0.1*10^4;
% thresh=2.2*10^4; % better used for vgat theta
% % thresh=0.05*10^4; % Used for theta
% thresh=0.015*10^4; % Used for vgat no theta
thresh=0.035*10^4; % Used for vgat theta
thresh=0.002*10^4;
thresh=0.01*10^4;
thresh=0.4*10^4;
thresh=0.44*10^4; % for 4 Hz
% thresh=0.10*10^4; % for 14 Hz
% thresh=028*10^4; % for 20 Hz
thresh=0.6*10^4; % for 20 Hz
% thresh=0.4*10^4; % for 20 Hz
% thresh=0.39*10^4; % for 20 Hz
% thresh=0.30*10^4; % for 20 Hz
% thresh=0.021*10^4;
% thresh=0.029*10^4;
% % thresh=0.05*10^4;
% thresh=3.6*10^4;
% % % thresh=0;
% isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh & nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh) & nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>0 & nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>0;
% isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
% isBigEnough=nanmean(data1(:,t>=5 & t<=6),2)>thresh;
% isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh & nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
data1=data1(isBigEnough,:);
data2=data2(isBigEnough,:);
disp(sum(isBigEnough==1));












% norm=0;
% normWindow=[1 3];
% % normWindow=[3.5 4];
% % normWindow=[4 6.5];
% baseSubtract=1;
% % baseWindow=[0 2.5];
% % baseWindow=[0 3];
% % baseWindow=[3.9 4];
% % baseWindow=[0 0.7];
% % baseWindow=[0 0.5];
% % baseWindow=[0 0.34];
% baseWindow=[0 0.51];
% % baseWindow=[0.51 1]; % for LOW FREQ
% % baseWindow=[0.9 1]; % for HIGH FREQ
% baseWindow=[0 0.7];
% % baseWindow=[0 0.51]; 
% % baseWindow=[3.7 3.9];
% % baseWindow=[3 3.7];
% % baseWindow=[0 3];
% % baseWindow=[2.25 3.5];
% % ds=[];
% % ds=2;
% % ds=3;
% % ds=20;
% ds=1;
% doMed=0;
% doFill=1;
% 
% if length(trialDuration)>1
%     t=trialDuration;
% else
% %     t=linspace(0,trialDuration,size(data1,2));
% %     t=linspace(1,trialDuration-1,size(data1,2));
%     t=linspace(0.5,trialDuration-0.5,size(data1,2));
% %     t=linspace(0.0825,trialDuration-0.0825,size(data1,2));
% %     t=linspace(0.165,trialDuration-0.165,size(data1,2));
% end
% 
% 
% % % Test whether response is sustained
% % useUnit=zeros(1,size(data2,1));
% % for i=1:size(data2,1)
% %     currBase=data2(i,t>=0 & t<=3);
% %     currResponse=data2(i,t>=5 & t<=6);
% %     p=ranksum(currBase,currResponse);
% %     if p<0.05 && nanmean(currBase)*1<nanmean(currResponse)
% %         useUnit(i)=1;
% %     end
% % end
% % data1=data1(useUnit==1,:);
% % data2=data2(useUnit==1,:);
% 
% 
% if baseSubtract==1
%     for i=1:size(data1,1)
%         base1=nanmean(data1(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
%         base2=nanmean(data2(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
%         data1(i,:)=data1(i,:)-base1;
%         data2(i,:)=data2(i,:)-base2;
%     end
% end
% 
% % % % thresh=0.15*10^4;
% % % % thresh=0.1*10^4;
% % % thresh=2.2*10^4; % better used for vgat theta
% % % % thresh=0.05*10^4; % Used for theta
% % % thresh=0.015*10^4; % Used for vgat no theta
% % thresh=0.035*10^4; % Used for vgat theta
% % thresh=0.002*10^4;
% % thresh=0.01*10^4;
% % thresh=0.4*10^4;
% % % thresh=0.44*10^4; % for 4 Hz
% % % thresh=0.35*10^4; % for 4 Hz
% % % thresh=0.26*10^4; % for 6 Hz
% % % thresh=0.13*10^4; % for 8 Hz
% % % thresh=0.27*10^4; % for 10 Hz
% % % thresh=0.4*10^4; % for 12 Hz
% % % thresh=0.15*10^4; % for 14 Hz
% % % thresh=0.10*10^4; % for 14 Hz
% % % thresh=0.6*10^4; % for 20 Hz
% % % thresh=0.1*10^4; % for 16 Hz
% % % thresh=0.2*10^4; % for 20 Hz
% % % thresh=0.7*10^4; % for 30 Hz
% % % thresh=0.7*10^4; % for 30 Hz
% % % thresh=0.9*10^4; % for 16 Hz
% % % thresh=1.2*10^4; % for 20 Hz theta vs no theta
% % thresh=0*10^4; % for TF tuning HIGH FREQ
% % % thresh=0.021*10^4;
% % % thresh=0.029*10^4;
% % % % thresh=0.05*10^4;
% % % thresh=3.6*10^4;
% % % % % thresh=0;
% % % isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh & nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh) & nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>0 & nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>0;
% % % isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
% % % isBigEnough=nanmean(data1(:,t>=5 & t<=6),2)>thresh;
% % % isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh & nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
% % % isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh | nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
% % % isBigEnough=(nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
% % % isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
% % data1=data1(isBigEnough,:);
% % data2=data2(isBigEnough,:);
% % disp(sum(isBigEnough==1));


if norm==1
    for i=1:size(data1,1)
        scale=nanmean(data1(i,t>=normWindow(1) & t<=normWindow(2)),2);
        data1(i,:)=data1(i,:)./scale;
%         scale=nanmean(data2(i,t>=normWindow(1) & t<=normWindow(2)),2);
        data2(i,:)=data2(i,:)./scale;
    end
end
%   
% if baseSubtract==1
%     for i=1:size(data1,1)
%         base1=nanmean(data1(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
%         base2=nanmean(data2(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
%         data1(i,:)=data1(i,:)-base1;
%         data2(i,:)=data2(i,:)-base2;
%     end
% end

% thresh=5*10^4;
% isBigEnough=nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh | nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh;
% data1=data1(isBigEnough,:);
% data2=data2(isBigEnough,:);

if ~isempty(ds)
    data1=downSampMatrix(data1,ds);
    data2=downSampMatrix(data2,ds);
    if length(trialDuration)>1
        trialDuration=downSampAv(trialDuration,ds);
    else
        trialDuration=downSampAv(t,ds);
    end
end

if doMed==1
    if length(trialDuration)>1
        figure();
        plot(trialDuration,nanmedian(data1,1),'Color',c1);
        hold on;
        if doFill==1
            fill([trialDuration fliplr(trialDuration)],[prctile(data1,55) fliplr(prctile(data1,45))],[0.5 0.5 0.5]);
        end
        plot(trialDuration,nanmedian(data1,1),'Color',c1);
%         plot(trialDuration,nanmedian(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
%         plot(trialDuration,nanmedian(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
        
        plot(trialDuration,nanmedian(data2,1),'Color',c2);
        if doFill==1
            fill([trialDuration fliplr(trialDuration)],[prctile(data2,55) fliplr(prctile(data2,45))],[0.1 0.7 0.5]);
        end
        plot(trialDuration,nanmedian(data2,1),'Color',c2);
%         plot(trialDuration,nanmedian(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
%         plot(trialDuration,nanmedian(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
        
        plot(trialDuration,nanmedian(data1,1),'Color',c1);
    else
        figure();
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1),'Color',c1);
        hold on;
        if doFill==1
            fill([linspace(0,trialDuration,size(data1,2)) fliplr(linspace(0,trialDuration,size(data1,2)))],[nanmedian(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)) fliplr(nanmedian(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)))],[0.5 0.5 0.5]);
        end
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1),'Color',c1);
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
        
        plot(linspace(0,trialDuration,size(data2,2)),nanmedian(data2,1),'Color',c2);
        if doFill==1
            fill([linspace(0,trialDuration,size(data2,2)) fliplr(linspace(0,trialDuration,size(data2,2)))],[nanmedian(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)) fliplr(nanmedian(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)))],[0.1 0.7 0.5]);
        end
        plot(linspace(0,trialDuration,size(data2,2)),nanmedian(data2,1),'Color',c2);
        plot(linspace(0,trialDuration,size(data2,2)),nanmedian(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
        plot(linspace(0,trialDuration,size(data2,2)),nanmedian(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
        
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1),'Color',c1);
    end
    return
end

if length(trialDuration)>1
    dropthebase=[0 1];
%     dropthebase=[3 3.7];
%     dropthebase=[0 4];
%     dropthebase=[0 3];
%     average_baseline_data1=nanmean(nanmean(data1(:,trialDuration>=dropthebase(1) & trialDuration<=dropthebase(2)),2),1);
%     average_baseline_data2=nanmean(nanmean(data2(:,trialDuration>=dropthebase(1) & trialDuration<=dropthebase(2)),2),1);
    average_baseline_data1=0;
    average_baseline_data2=0;
    figure(); 
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1,'Color',c1);
    hold on;
    if doFill==1
        fill([trialDuration fliplr(trialDuration)],[nanmean(data1,1)-average_baseline_data1+nanstd(data1,[],1)./sqrt(size(data1,1)) fliplr(nanmean(data1,1)-average_baseline_data1-nanstd(data1,[],1)./sqrt(size(data1,1)))],[0.5 0.5 0.5]);
    end
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1,'Color',c1);
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1+nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1-nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
    
    plot(trialDuration,nanmean(data2,1)-average_baseline_data2,'Color',c2);
    if doFill==1
        fill([trialDuration fliplr(trialDuration)],[nanmean(data2,1)-average_baseline_data2+nanstd(data2,[],1)./sqrt(size(data2,1)) fliplr(nanmean(data2,1)-average_baseline_data2-nanstd(data2,[],1)./sqrt(size(data2,1)))],[0.1 0.7 0.5]);
    end
    plot(trialDuration,nanmean(data2,1)-average_baseline_data2,'Color',c2);
    plot(trialDuration,nanmean(data2,1)-average_baseline_data2+nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
    plot(trialDuration,nanmean(data2,1)-average_baseline_data2-nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
    
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1,'Color',c1);
    
    disp('Evoked no LED no theta');
    disp(nanmean(nanmean(data1(:,trialDuration>=4 & trialDuration<=6.5),1)-average_baseline_data1));
    disp('Evoked LED no theta');
    disp(nanmean(nanmean(data2(:,trialDuration>=4 & trialDuration<=6.5),1)-average_baseline_data2));
    
%     figure(); 
%     plot(trialDuration,data1','Color',c1);
%     hold on;
%     plot(trialDuration,data2','Color',c2);
    
else
    figure();
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1),'Color',c1);
    hold on;
    if doFill==1
        fill([linspace(0,trialDuration,size(data1,2)) fliplr(linspace(0,trialDuration,size(data1,2)))],[nanmean(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)) fliplr(nanmean(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)))],[0.5 0.5 0.5]);
    end
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1),'Color',c1);
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
    
    plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1),'Color',c2);
    if doFill==1
        fill([linspace(0,trialDuration,size(data2,2)) fliplr(linspace(0,trialDuration,size(data2,2)))],[nanmean(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)) fliplr(nanmean(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)))],[0.1 0.7 0.5]);
    end
    plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1),'Color',c2);
    plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
    plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
    
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1),'Color',c1);
    
end
end