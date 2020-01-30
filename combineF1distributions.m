function combineF1distributions(datadir,addName,suppressOutput)

responsiveOnly=true;
visResponseThresh_backup='50';
onlyTakeAboveZero=true;

% bincenters=[-10:0.01:10];
bincenters=[-100:0.1:100];
nbins=bincenters;

warning('off')

nconditions=1;
if ~isempty(addName)
    if ischar(addName)
        for i=1:length(datadir)
            temp=datadir{i};
            datadir{i}=[temp addName];
        end
        nconditions=1;
    elseif iscell(addName)
        for j=1:length(addName)
            currAddName=addName{j};
            for i=1:length(datadir)
                temp=datadir{i};
                newdatadir{length(datadir)*(j-1)+i}=[temp currAddName];
            end
        end
        datadir=newdatadir;
        nconditions=length(addName);
    end
end

if iscell(datadir)
    
    all_noTheta_noLED_F1=[];
    all_noTheta_LED_F1=[];
    all_theta_noLED_F1=[];
    all_theta_LED_F1=[];
    
    all_noTheta_noLED_F1_zscoreNow=[];
    all_noTheta_LED_F1_zscoreNow=[];
    all_theta_noLED_F1_zscoreNow=[];
    all_theta_LED_F1_zscoreNow=[];
    
    all_noTheta_noLED_alphas=[];
    all_noTheta_LED_alphas=[];
    all_theta_noLED_alphas=[];
    all_theta_LED_alphas=[];
    
    all_noTheta_noLED_F1_hists=[];
    all_noTheta_LED_F1_hists=[];
    all_theta_noLED_F1_hists=[];
    all_theta_LED_F1_hists=[];
    
    all_noTheta_noLED_alphas_hists=[];
    all_noTheta_LED_alphas_hists=[];
    all_theta_noLED_alphas_hists=[];
    all_theta_LED_alphas_hists=[];
    
    all_noTheta_noLED_F1_forced=[];
    all_noTheta_LED_F1_forced=[];
    all_theta_noLED_F1_forced=[];
    all_theta_LED_F1_forced=[];
    
    all_noTheta_noLED_alphas_forced=[];
    all_noTheta_LED_alphas_forced=[];
    all_theta_noLED_alphas_forced=[];
    all_theta_LED_alphas_forced=[];
    
    nused=0;
    
    for i=1:length(datadir)
        visResponseThresh=visResponseThresh_backup;
        
        d=datadir{i};
        
        if ~exist([d '\noTheta_noLED_out.mat'],'file')
            % for DG only
            [noTheta_noLED_out,noTheta_LED_out,theta_noLED_out,theta_LED_out]=getF1distributions(d,[],[],[],[],[2.25 3.75],{[10 14],[1 50]},[4.5 6],[3.7 3.9],true,[],[],[],[]);
            save([d '\noTheta_noLED_out.mat']);
            save([d '\noTheta_LED_out.mat']);
            save([d '\theta_noLED_out.mat']);
            save([d '\theta_LED_out.mat']);
        else
            a=load([d '\' 'noTheta_noLED_out.mat']);
            noTheta_noLED_out=a.noTheta_noLED_out;
            a=load([d '\' 'noTheta_LED_out.mat']);
            noTheta_LED_out=a.noTheta_LED_out;
            a=load([d '\' 'theta_noLED_out.mat']);
            theta_noLED_out=a.theta_noLED_out;
            a=load([d '\' 'theta_LED_out.mat']);
            theta_LED_out=a.theta_LED_out;
        end
        
        if responsiveOnly==true
             visResponses=nanmean(noTheta_noLED_out.F1s,2) + nanmean(noTheta_LED_out.F1s,2) + nanmean(theta_noLED_out.F1s,2) + nanmean(theta_LED_out.F1s,2);
             if ischar(visResponseThresh)
                 visResponseThresh=prctile(visResponses,str2num(visResponseThresh));
                 nused=nused+nansum(visResponses>visResponseThresh);
             end
             noTheta_noLED_out.F1s=noTheta_noLED_out.F1s(visResponses>visResponseThresh,:);
             noTheta_LED_out.F1s=noTheta_LED_out.F1s(visResponses>visResponseThresh,:);
             theta_noLED_out.F1s=theta_noLED_out.F1s(visResponses>visResponseThresh,:);
             theta_LED_out.F1s=theta_LED_out.F1s(visResponses>visResponseThresh,:);
             
             noTheta_noLED_out.alphas=noTheta_noLED_out.alphas(visResponses>visResponseThresh,:);
             noTheta_LED_out.alphas=noTheta_LED_out.alphas(visResponses>visResponseThresh,:);
             theta_noLED_out.alphas=theta_noLED_out.alphas(visResponses>visResponseThresh,:);
             theta_LED_out.alphas=theta_LED_out.alphas(visResponses>visResponseThresh,:);
             
             if onlyTakeAboveZero==true
                 noTheta_noLED_out.F1s(noTheta_noLED_out.F1s<0)=nan;
                 noTheta_LED_out.F1s(noTheta_LED_out.F1s<0)=nan;
                 theta_noLED_out.F1s(theta_noLED_out.F1s<0)=nan;
                 theta_LED_out.F1s(theta_LED_out.F1s<0)=nan;
             end
        end
        
        z_noTheta_noLED_out=zscore_data(noTheta_noLED_out.F1s,[noTheta_noLED_out.F1s noTheta_LED_out.F1s theta_noLED_out.F1s theta_LED_out.F1s]);
        z_noTheta_LED_out=zscore_data(noTheta_LED_out.F1s,[noTheta_noLED_out.F1s noTheta_LED_out.F1s theta_noLED_out.F1s theta_LED_out.F1s]);
        z_theta_noLED_out=zscore_data(theta_noLED_out.F1s,[noTheta_noLED_out.F1s noTheta_LED_out.F1s theta_noLED_out.F1s theta_LED_out.F1s]);
        z_theta_LED_out=zscore_data(theta_LED_out.F1s,[noTheta_noLED_out.F1s noTheta_LED_out.F1s theta_noLED_out.F1s theta_LED_out.F1s]);
        
        all_noTheta_noLED_F1_zscoreNow=[all_noTheta_noLED_F1_zscoreNow reshape(z_noTheta_noLED_out(1:end),1,length(z_noTheta_noLED_out(1:end)))];
        all_noTheta_LED_F1_zscoreNow=[all_noTheta_LED_F1_zscoreNow reshape(z_noTheta_LED_out(1:end),1,length(z_noTheta_LED_out(1:end)))];
        all_theta_noLED_F1_zscoreNow=[all_theta_noLED_F1_zscoreNow reshape(z_theta_noLED_out(1:end),1,length(z_theta_noLED_out(1:end)))];
        all_theta_LED_F1_zscoreNow=[all_theta_LED_F1_zscoreNow reshape(z_theta_LED_out(1:end),1,length(z_theta_LED_out(1:end)))];

        all_noTheta_noLED_F1=[all_noTheta_noLED_F1 reshape(noTheta_noLED_out.F1s(1:end),1,length(noTheta_noLED_out.F1s(1:end)))];
        all_noTheta_LED_F1=[all_noTheta_LED_F1 reshape(noTheta_LED_out.F1s(1:end),1,length(noTheta_LED_out.F1s(1:end)))];
        all_theta_noLED_F1=[all_theta_noLED_F1 reshape(theta_noLED_out.F1s(1:end),1,length(theta_noLED_out.F1s(1:end)))];
        all_theta_LED_F1=[all_theta_LED_F1 reshape(theta_LED_out.F1s(1:end),1,length(theta_LED_out.F1s(1:end)))];
        
        all_noTheta_noLED_alphas=[all_noTheta_noLED_alphas reshape(noTheta_noLED_out.alphas(1:end),1,length(noTheta_noLED_out.alphas(1:end)))];
        all_noTheta_LED_alphas=[all_noTheta_LED_alphas reshape(noTheta_LED_out.alphas(1:end),1,length(noTheta_LED_out.alphas(1:end)))];
        all_theta_noLED_alphas=[all_theta_noLED_alphas reshape(theta_noLED_out.alphas(1:end),1,length(theta_noLED_out.alphas(1:end)))];
        all_theta_LED_alphas=[all_theta_LED_alphas reshape(theta_LED_out.alphas(1:end),1,length(theta_LED_out.alphas(1:end)))];
        
        for j=1:size(noTheta_noLED_out.F1s,1)
            all_noTheta_noLED_F1_hists=[all_noTheta_noLED_F1_hists; hist(noTheta_noLED_out.F1s(j,:),bincenters)];
            all_noTheta_LED_F1_hists=[all_noTheta_LED_F1_hists; hist(noTheta_LED_out.F1s(j,:),bincenters)];
            all_theta_noLED_F1_hists=[all_theta_noLED_F1_hists; hist(theta_noLED_out.F1s(j,:),bincenters)];
            all_theta_LED_F1_hists=[all_theta_LED_F1_hists; hist(theta_LED_out.F1s(j,:),bincenters)];
            
            all_noTheta_noLED_alphas_hists=[all_noTheta_noLED_alphas_hists; hist(noTheta_noLED_out.alphas(j,:),bincenters)];
            all_noTheta_LED_alphas_hists=[all_noTheta_LED_alphas_hists; hist(noTheta_LED_out.alphas(j,:),bincenters)];
            all_theta_noLED_alphas_hists=[all_theta_noLED_alphas_hists; hist(theta_noLED_out.alphas(j,:),bincenters)];
            all_theta_LED_alphas_hists=[all_theta_LED_alphas_hists; hist(theta_LED_out.alphas(j,:),bincenters)];
            
            all_noTheta_noLED_F1_forced=[all_noTheta_noLED_F1_forced forceTo2Gaussians(noTheta_noLED_out.F1s(j,:))];
            all_noTheta_LED_F1_forced=[all_noTheta_LED_F1_forced forceTo2Gaussians(noTheta_LED_out.F1s(j,:))];
            all_theta_noLED_F1_forced=[all_theta_noLED_F1_forced forceTo2Gaussians(theta_noLED_out.F1s(j,:))];
            all_theta_LED_F1_forced=[all_theta_LED_F1_forced forceTo2Gaussians(theta_LED_out.F1s(j,:))];
            
            all_noTheta_noLED_alphas_forced=[all_noTheta_noLED_alphas_forced forceTo2Gaussians(noTheta_noLED_out.alphas(j,:))];
            all_noTheta_LED_alphas_forced=[all_noTheta_LED_alphas_forced forceTo2Gaussians(noTheta_LED_out.alphas(j,:))];
            all_theta_noLED_alphas_forced=[all_theta_noLED_alphas_forced forceTo2Gaussians(theta_noLED_out.alphas(j,:))];
            all_theta_LED_alphas_forced=[all_theta_LED_alphas_forced forceTo2Gaussians(theta_LED_out.alphas(j,:))];
        end
      
    end
else
    error('expected datadir to be cell array');
end

if suppressOutput==false
    disp('n');
    disp(nused);
    disp('n conditions');
    disp(nconditions);
    
    [outx_con,outn_con]=plotDistributions(all_noTheta_noLED_F1,all_noTheta_noLED_alphas,'no theta no LED',downSampAv(nbins,4));
    [outx_led,outn_led]=plotDistributions(all_noTheta_LED_F1,all_noTheta_LED_alphas,'no theta LED',downSampAv(nbins,4));
    disp('nanmean no theta');
    disp(nanmean(all_noTheta_noLED_F1(all_noTheta_noLED_F1>0)));
    disp('nanmean no theta LED');
    disp(nanmean(all_noTheta_LED_F1(all_noTheta_LED_F1>0)));
%     plotDistributions(all_theta_noLED_F1,all_theta_noLED_alphas,'theta no LED',downSampAv(nbins,8));
%     plotDistributions(all_theta_LED_F1,all_theta_LED_alphas,'theta LED',downSampAv(nbins,8));  
    figure(); plot(outx_con,outn_con./nansum(outn_con),'Color','k'); hold on; plot(outx_led,outn_led./nansum(outn_led),'Color','b'); title('no theta con vs led');
    
    [outx_con,outn_con]=plotDistributions(all_noTheta_noLED_F1_zscoreNow,all_noTheta_noLED_alphas,'no theta no LED zscore',downSampAv([-10:0.05:10],1));
    [outx_led,outn_led]=plotDistributions(all_noTheta_LED_F1_zscoreNow,all_noTheta_LED_alphas,'no theta LED zscore',downSampAv([-10:0.05:10],1));
    figure(); plot(outx_con,outn_con./nansum(outn_con),'Color','k'); hold on; plot(outx_led,outn_led./nansum(outn_led),'Color','b'); title('no theta con vs led zscore');
    figure(); plot(outx_con,smooth(outn_con./nansum(outn_con)-outn_led./nansum(outn_led),50),'Color','r'); title('control minus led');

    
    [outx_con,outn_con]=plotDistributions(all_theta_noLED_F1_zscoreNow,all_theta_noLED_alphas,'theta no LED zscore',downSampAv([-10:0.05:10],2));
    [outx_led,outn_led]=plotDistributions(all_theta_LED_F1_zscoreNow,all_theta_LED_alphas,'theta LED zscore',downSampAv([-10:0.05:10],2));  
    figure(); plot(outx_con,outn_con./nansum(outn_con),'Color','k'); hold on; plot(outx_led,outn_led./nansum(outn_led),'Color','b'); title('theta con vs led zscore');
    figure(); plot(outx_con,smooth(outn_con./nansum(outn_con)-outn_led./nansum(outn_led),50),'Color','r'); title('control minus led theta');

%     plotDistributions(all_noTheta_noLED_F1_forced,all_noTheta_noLED_alphas_forced,'no theta no LED forced to 2 gauss',downSampAv(nbins,5));
%     plotDistributions(all_noTheta_LED_F1_forced,all_noTheta_LED_alphas_forced,'no theta LED forced to 2 gauss',downSampAv(nbins,5));
%     plotDistributions(all_theta_noLED_F1_forced,all_theta_noLED_alphas_forced,'theta no LED forced to 2 gauss',downSampAv(nbins,5));
%     plotDistributions(all_theta_LED_F1_forced,all_theta_LED_alphas_forced,'theta LED forced to 2 gauss',downSampAv(nbins,5));    
end

end

function F1s_zscored=zscore_data(F1s_per_unit,popData)

% Z-score (sample - mean)/stdev
F1s_zscored=nan(size(F1s_per_unit));
for i=1:size(F1s_per_unit,1)
    F1s_zscored(i,:)=(F1s_per_unit(i,:)-nanmean(popData(i,:),2))./nanstd(popData(i,:),[],2);
end

end

function out_data=forceTo2Gaussians(data)

if nansum(~isnan(data))<=1
    out_data=nan(size(data));
    return
end

try
    gmmodel=fitgmdist(data',2,'RegularizationValue',0.01);
    mus=gmmodel.mu;
    mus=sort(mus);
    diffBetweenMus=mus(2)-mus(1);
    scaleFac=2/diffBetweenMus;
    out_data=data*scaleFac;
    bottomMu=mus(1)*scaleFac;
    out_data=out_data-bottomMu-1;
catch
    out_data=nan(size(data));
end

end

function [outx,outn]=plotDistributions(temp,temp_alpha,tit,nbins)

[n,x]=hist(temp,nbins);
outx=x;
outn=n;
figure();
plot(x,n,'Color','k');
maxn=nanmax(n);
hold on;
if length(nbins)==1
    smallerbins=floor(nbins/4);
else
    smallerbins=downSampAv(nbins,4);
end
[n,x]=hist(temp(temp_alpha<1),smallerbins);
plot(x,(n./nanmax(n))*maxn,'Color','g');
[n,x]=hist(temp(temp_alpha>1),smallerbins);
plot(x,(n./nanmax(n))*maxn,'Color','r');
title(tit);

end