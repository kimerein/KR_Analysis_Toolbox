function combineSpatFreqs(datadir)

allSpatFreqs=[0.02 0.03 0.04 0.06 0.08];

all_noTheta_noLED_numer=nan(1000,length(allSpatFreqs));
all_noTheta_noLED_denom=nan(1000,length(allSpatFreqs));
all_theta_noLED_numer=nan(1000,length(allSpatFreqs));
all_theta_noLED_denom=nan(1000,length(allSpatFreqs));
all_noTheta_LED_numer=nan(1000,length(allSpatFreqs));
all_theta_LED_numer=nan(1000,length(allSpatFreqs));

ls=dir(datadir);
k=1;
for i=3:length(ls)
    disp(i);
    % read in each day's data
    a=load([ls(i).folder '\' ls(i).name '\spatFreqs.mat']);
    spatFreqs=a.spatFreqs;
    % find indices into allSpatFreqs
    indIntoAll=nan(1,length(spatFreqs));
    for j=1:length(spatFreqs)
        [~,indIntoAll(j)]=nanmin(abs(spatFreqs(j)-allSpatFreqs));
    end
    if exist([ls(i).folder '\' ls(i).name '\useTheseUnits.mat'],'file')
        a=load([ls(i).folder '\' ls(i).name '\useTheseUnits.mat']);
        useTheseUnits=a.useTheseUnits;
    else
        useTheseUnits=[];
    end
    a=load([ls(i).folder '\' ls(i).name '\noTheta_noLED.mat']);
    if ~isempty(useTheseUnits)
        a.numer=a.numer(useTheseUnits==0,:);
        a.denom=a.denom(useTheseUnits==0,:);
    end
    for j=1:length(indIntoAll)
        all_noTheta_noLED_numer(k:k+(size(a.numer,1)-1),indIntoAll(j))=a.numer(:,j);
        all_noTheta_noLED_denom(k:k+(size(a.numer,1)-1),indIntoAll(j))=a.denom(:,j);
    end
    a=load([ls(i).folder '\' ls(i).name '\theta_noLED.mat']);
    if ~isempty(useTheseUnits)
        a.numer=a.numer(useTheseUnits==0,:);
        a.denom=a.denom(useTheseUnits==0,:);
    end
    for j=1:length(indIntoAll)
        all_theta_noLED_numer(k:k+(size(a.numer,1)-1),indIntoAll(j))=a.numer(:,j);
        all_theta_noLED_denom(k:k+(size(a.numer,1)-1),indIntoAll(j))=a.denom(:,j);
    end
    a=load([ls(i).folder '\' ls(i).name '\noTheta_LED.mat']);
    if ~isempty(useTheseUnits)
        a.numer=a.numer(useTheseUnits==0,:);
        a.denom=a.denom(useTheseUnits==0,:);
    end
    for j=1:length(indIntoAll)
        all_noTheta_LED_numer(k:k+(size(a.numer,1)-1),indIntoAll(j))=a.numer(:,j);
    end
    a=load([ls(i).folder '\' ls(i).name '\theta_LED.mat']);
    if ~isempty(useTheseUnits)
        a.numer=a.numer(useTheseUnits==0,:);
        a.denom=a.denom(useTheseUnits==0,:);
    end
    for j=1:length(indIntoAll)
        all_theta_LED_numer(k:k+(size(a.numer,1)-1),indIntoAll(j))=a.numer(:,j);
    end
    k=k+size(a.numer,1);
end

sz=100;

all_theta_LED_numer(~isnan(all_noTheta_noLED_numer(:,4)),:)=nan(nansum(~isnan(all_noTheta_noLED_numer(:,4))),5);
all_theta_noLED_numer(~isnan(all_noTheta_noLED_numer(:,4)),:)=nan(nansum(~isnan(all_noTheta_noLED_numer(:,4))),5);
all_noTheta_LED_numer(~isnan(all_noTheta_noLED_numer(:,4)),:)=nan(nansum(~isnan(all_noTheta_noLED_numer(:,4))),5);
all_noTheta_noLED_numer(~isnan(all_noTheta_noLED_numer(:,4)),:)=nan(nansum(~isnan(all_noTheta_noLED_numer(:,4))),5);

% discardNonResponsive=all_noTheta_noLED_numer(1:end)<0 | all_theta_noLED_numer(1:end)<0;
% all_theta_LED_numer(discardNonResponsive)=nan; 
% all_noTheta_LED_numer(discardNonResponsive)=nan; 
% all_theta_noLED_numer(discardNonResponsive)=nan; 
% all_noTheta_noLED_numer(discardNonResponsive)=nan; 
% disp(nansum(~discardNonResponsive));

currsum=0;
for i=[1 2 3 5]
    tempdata=all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
    currsum=currsum+nanmean(tempdata(~isinf(tempdata)));
end
currsum=currsum./length(allSpatFreqs);

i=1;
figure(); 
r = 0 + (0.2-0) .* rand(1000,1);
% scatter(ones(1000,1)+r,all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(1.2*ones(1000,1)+r,all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
line([1-0.025 1+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([1 1],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
line([1.05-0.025 1.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([1.05 1.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);


i=2;
% scatter(2*ones(1000,1)+r,all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(2.2*ones(1000,1)+r,all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
line([2-0.025 2+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([2 2],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
line([2.05-0.025 2.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([2.05 2.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=3;
% scatter(3*ones(1000,1)+r,all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(3.2*ones(1000,1)+r,all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
line([3-0.025 3+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([3 3],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
tempdata=tempdata(~isinf(tempdata));
line([3.05-0.025 3.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([3.05 3.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);


i=4;
% scatter(4*ones(1000,1)+r,all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(4.2*ones(1000,1)+r,all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
line([4-0.025 4+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([4 4],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
line([4.05-0.025 4.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([4.05 4.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=5;
% scatter(5*ones(1000,1)+r,all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(5.2*ones(1000,1)+r,all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
line([5-0.025 5+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([5 5],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2))./currsum;
line([5.05-0.025 5.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([5.05 5.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);
title('No theta; norm to con all spat freqs');


% Now theta

currsum=0;
for i=[1 2 3 5]
    tempdata=all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2);
    currsum=currsum+nanmean(tempdata(~isinf(tempdata)));
end
currsum=currsum./length(allSpatFreqs);

i=1;
figure(); 
% scatter(ones(1000,1)+r,all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(1.2*ones(1000,1)+r,all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
line([1-0.025 1+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([1 1],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
line([1.05-0.025 1.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([1.05 1.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=2;
% scatter(2*ones(1000,1)+r,all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(2.2*ones(1000,1)+r,all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
line([2-0.025 2+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([2 2],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
line([2.05-0.025 2.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([2.05 2.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=3;
% scatter(3*ones(1000,1)+r,all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(3.2*ones(1000,1)+r,all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
line([3-0.025 3+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([3 3],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
tempdata=tempdata(~isinf(tempdata));
line([3.05-0.025 3.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([3.05 3.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=4;
% scatter(4*ones(1000,1)+r,all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(4.2*ones(1000,1)+r,all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
line([4-0.025 4+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([4 4],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
line([4.05-0.025 4.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([4.05 4.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=5;
% scatter(5*ones(1000,1)+r,all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(5.2*ones(1000,1)+r,all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=(all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
tempdata=tempdata(~isinf(tempdata));
line([5-0.025 5+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([5 5],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=(all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2))./currsum;
line([5.05-0.025 5.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([5.05 5.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);
title('Theta; norm to con all spat freqs');

figure();
scatter(all_theta_noLED_numer(:,i)./nanmean(all_theta_noLED_denom,2),all_theta_LED_numer(:,i)./nanmean(all_theta_noLED_denom,2));


% No theta vs theta
% discardNonResponsive=all_noTheta_noLED_numer(1:end)<0 | all_theta_noLED_numer(1:end)<0;
% all_theta_LED_numer(discardNonResponsive)=nan; 
% all_noTheta_LED_numer(discardNonResponsive)=nan; 
% all_theta_noLED_numer(discardNonResponsive)=nan; 
% all_noTheta_noLED_numer(discardNonResponsive)=nan; 


i=1;
figure(); 
scatter(ones(1000,1),all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
hold on; scatter(1.05*ones(1000,1),all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([1-0.025 1+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([1 1],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([1.05-0.025 1.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([1.05 1.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=2;
scatter(2*ones(1000,1),all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
hold on; scatter(2.05*ones(1000,1),all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([2-0.025 2+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([2 2],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([2.05-0.025 2.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([2.05 2.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=3;
scatter(3*ones(1000,1),all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
hold on; scatter(3.05*ones(1000,1),all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([3-0.025 3+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([3 3],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
tempdata=tempdata(~isinf(tempdata));
line([3.05-0.025 3.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([3.05 3.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=4;
scatter(4*ones(1000,1),all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
hold on; scatter(4.05*ones(1000,1),all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([4-0.025 4+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([4 4],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([4.05-0.025 4.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([4.05 4.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=5;
scatter(5*ones(1000,1),all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
hold on; scatter(5.05*ones(1000,1),all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([5-0.025 5+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([5 5],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_noLED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([5.05-0.025 5.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([5.05 5.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);
title('No theta vs theta; norm to no theta all spat freqs');

% No theta LED vs theta LED
i=1;
figure(); 
% scatter(ones(1000,1),all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(1.05*ones(1000,1),all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([1-0.025 1+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([1 1],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([1.05-0.025 1.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([1.05 1.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=2;
% scatter(2*ones(1000,1),all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(2.05*ones(1000,1),all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([2-0.025 2+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([2 2],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([2.05-0.025 2.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([2.05 2.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=3;
% scatter(3*ones(1000,1),all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(3.05*ones(1000,1),all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([3-0.025 3+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([3 3],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
tempdata=tempdata(~isinf(tempdata));
line([3.05-0.025 3.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([3.05 3.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=4;
% scatter(4*ones(1000,1),all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(4.05*ones(1000,1),all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([4-0.025 4+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([4 4],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([4.05-0.025 4.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([4.05 4.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);

i=5;
% scatter(5*ones(1000,1),all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,'k'); %alpha 0.05;
% hold on; scatter(5.05*ones(1000,1),all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2),sz,[0    0.5938    1.0000]); %alpha 0.05;
tempdata=all_noTheta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([5-0.025 5+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color','k');
line([5 5],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color','k');
tempdata=all_theta_LED_numer(:,i)./nanmean(all_noTheta_noLED_denom,2);
line([5.05-0.025 5.05+0.025],[nanmean(tempdata) nanmean(tempdata)],'Color',[0    0.5938    1.0000]);
line([5.05 5.05],[nanmean(tempdata)-nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata))) nanmean(tempdata)+nanstd(tempdata,[],1)./sqrt(nansum(~isnan(tempdata)))],'Color',[0    0.5938    1.0000]);
title('No theta LED vs theta LED; norm to no theta con all spat freqs');


end