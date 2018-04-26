function putTogetherGamma(useDir,listing)
    
allGamma.autocorr_x=[];
allGamma.autocorr_y=[];
allGamma.autocorr_x_led=[];
allGamma.autocorr_y_led=[];
allGamma.power_x=[];
allGamma.power_y=[];
allGamma.power_x_led=[];
allGamma.power_y_led=[];

allNo.autocorr_x=[];
allNo.autocorr_y=[];
allNo.autocorr_x_led=[];
allNo.autocorr_y_led=[];
allNo.power_x=[];
allNo.power_y=[];
allNo.power_x_led=[];
allNo.power_y_led=[];
for i=3:length(listing)
    currFolder=listing(i).name;
    l=dir([useDir '\' currFolder '\']);
    for j=3:length(l)
        a.gamma=[];
        a.no_gamma=[];
        a=load([useDir '\' currFolder '\' l(j).name]);
        if isfield(a,'gamma')
            gamma=a.gamma;
            no_gamma=a.no_gamma;
        else
%             gamma=a.allUnitGamma;
            gamma=a.allTogether2;
            no_gamma=[];
        end
        for k=1:length(gamma)
            g=gamma(k);
            allGamma.autocorr_x=[allGamma.autocorr_x; g.autocorr_x];
            allGamma.autocorr_y=[allGamma.autocorr_y; g.autocorr_y];
            allGamma.autocorr_x_led=[allGamma.autocorr_x_led; g.autocorr_x_led];
            allGamma.autocorr_y_led=[allGamma.autocorr_y_led; g.autocorr_y_led];
            allGamma.power_x=[allGamma.power_x; g.power_x];
            allGamma.power_y=[allGamma.power_y; g.power_y'];
            allGamma.power_x_led=[allGamma.power_x_led; g.power_x_led];
            allGamma.power_y_led=[allGamma.power_y_led; g.power_y_led'];
        end
        for k=1:length(no_gamma)
            g=no_gamma(k);
            allNo.autocorr_x=[allNo.autocorr_x; g.autocorr_x];
            allNo.autocorr_y=[allNo.autocorr_y; g.autocorr_y];
            allNo.autocorr_x_led=[allNo.autocorr_x_led; g.autocorr_x_led];
            allNo.autocorr_y_led=[allNo.autocorr_y_led; g.autocorr_y_led];
            allNo.power_x=[allNo.power_x; g.power_x];
            allNo.power_y=[allNo.power_y; g.power_y'];
            allNo.power_x_led=[allNo.power_x_led; g.power_x_led];
            allNo.power_y_led=[allNo.power_y_led; g.power_y_led'];
        end
    end
end

% Normalize then plot all power spec of auto-corr
temp=allGamma.power_y;
% allGamma.power_y=allGamma.power_y./repmat(max(allGamma.power_y,[],2),1,size(allGamma.power_y,2));
% allGamma.power_y_led=allGamma.power_y_led./repmat(max(temp,[],2),1,size(temp,2));

figure(); 
plot(nanmean(allGamma.power_x,1),allGamma.power_y','Color','k');
figure();
plot(nanmean(allGamma.power_x_led,1),allGamma.power_y_led','Color','r');

figure(); 
hax=axes();
hl=plotLineAndErr(nanmean(allGamma.power_x,1),allGamma.power_y,'k',hax);
hold on;
hl=plotLineAndErr(nanmean(allGamma.power_x_led,1),allGamma.power_y_led,'r',hax);
title('Gamma Units Power Spec');

temp=allNo.power_y;
% allNo.power_y=allNo.power_y./repmat(max(allNo.power_y,[],2),1,size(allNo.power_y,2));
% allNo.power_y_led=allNo.power_y_led./repmat(max(temp,[],2),1,size(temp,2));

figure(); 
plot(nanmean(allNo.power_x,1),allNo.power_y','Color','k');
figure();
plot(nanmean(allNo.power_x_led,1),allNo.power_y_led','Color','r');

figure(); 
hax=axes();
hl=plotLineAndErr(nanmean(allNo.power_x,1),allNo.power_y,'k',hax);
hold on;
hl=plotLineAndErr(nanmean(allNo.power_x_led,1),allNo.power_y_led,'r',hax);
title('No Gamma Units Power Spec');

% Normalize then plot all auto-corr
temp=allGamma.autocorr_y;
allGamma.autocorr_y=allGamma.autocorr_y./repmat(max(allGamma.autocorr_y,[],2),1,size(allGamma.autocorr_y,2));
allGamma.autocorr_y_led=allGamma.autocorr_y_led./repmat(max(temp,[],2),1,size(temp,2));

figure(); 
plot(nanmean(allGamma.autocorr_x,1),allGamma.autocorr_y','Color','k');
figure();
plot(nanmean(allGamma.autocorr_x_led,1),allGamma.autocorr_y_led','Color','r');

figure(); 
hax=axes();
hl=plotLineAndErr(nanmean(allGamma.autocorr_x,1),allGamma.autocorr_y,'k',hax);
hold on;
hl=plotLineAndErr(nanmean(allGamma.autocorr_x_led,1),allGamma.autocorr_y_led,'r',hax);
title('Gamma Units Autocorr');

temp=allNo.autocorr_y;
allNo.autocorr_y=allNo.autocorr_y./repmat(max(allNo.autocorr_y,[],2),1,size(allNo.autocorr_y,2));
allNo.autocorr_y_led=allNo.autocorr_y_led./repmat(max(temp,[],2),1,size(temp,2));

figure(); 
plot(nanmean(allNo.autocorr_x,1),allNo.autocorr_y','Color','k');
figure();
plot(nanmean(allNo.autocorr_x_led,1),allNo.autocorr_y_led','Color','r');

figure(); 
hax=axes();
hl=plotLineAndErr(nanmean(allNo.autocorr_x,1),allNo.autocorr_y,'k',hax);
hold on;
hl=plotLineAndErr(nanmean(allNo.autocorr_x_led,1),allNo.autocorr_y_led,'r',hax);
title('No Gamma Units Autocorr');

% Normalize then plot all auto-corr (LED+no LED)
together.autocorr_y=[allGamma.autocorr_y; allNo.autocorr_y];
together.autocorr_y_led=[allGamma.autocorr_y_led; allNo.autocorr_y_led];
temp=together.autocorr_y;
together.autocorr_y=together.autocorr_y./repmat(max(together.autocorr_y,[],2),1,size(together.autocorr_y,2));
together.autocorr_y_led=together.autocorr_y_led./repmat(max(temp,[],2),1,size(temp,2));

figure(); 
plot(nanmean(allGamma.autocorr_x,1),together.autocorr_y','Color','k');
figure();
plot(nanmean(allGamma.autocorr_x_led,1),together.autocorr_y_led','Color','r');

figure(); 
hax=axes();
hl=plotLineAndErr(nanmean(allGamma.autocorr_x,1),together.autocorr_y,'k',hax);
hold on;
hl=plotLineAndErr(nanmean(allGamma.autocorr_x_led,1),together.autocorr_y_led,'r',hax);
title('Gamma and No Gamma Units Autocorr');


% Normalize then plot all power spec of auto-corr
together.power_y=[allGamma.power_y; allNo.power_y];
together.power_y_led=[allGamma.power_y_led; allNo.power_y_led];
temp=together.power_y;
% together.power_y=together.power_y./repmat(max(together.power_y,[],2),1,size(together.power_y,2));
% together.power_y_led=together.power_y_led./repmat(max(temp,[],2),1,size(temp,2));

figure(); 
plot(nanmean(allGamma.power_x,1),together.power_y','Color','k');
figure();
plot(nanmean(allGamma.power_x_led,1),together.power_y_led','Color','r');

figure(); 
hax=axes();
hl=plotLineAndErr(nanmean(allGamma.power_x,1),together.power_y,'k',hax);
hold on;
hl=plotLineAndErr(nanmean(allGamma.power_x_led,1),together.power_y_led,'r',hax);
title('Gamma and No Gamma Units Power Spec');

nbg=nanmean(together.power_y(:,allGamma.power_x(1,:)>=50 & allGamma.power_x(1,:)<=70),2);
nbg_led=nanmean(together.power_y_led(:,allGamma.power_x_led(1,:)>=50 & allGamma.power_x_led(1,:)<=70),2);
figure(); 
scatter(nbg,nbg_led);
xlabel('50-70 Hz Power');
ylabel('50-70 Hz Power With LED');

end

function hl=plotLineAndStd(x,data,c,hax)

hl=plot(x,nanmean(data,1),'Color',c);
addErrBar(x,nanmean(data,1),nanstd(data,[],1),'y',hax,hl);

end

function hl=plotLineAndErr(x,data,c,hax)

hl=plot(x,nanmean(data,1),'Color',c);
addErrBar(x,nanmean(data,1),nanstd(data,[],1)./sqrt(size(data,1)),'y',hax,hl);

end