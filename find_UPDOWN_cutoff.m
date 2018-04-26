function [fitted,window,led,histdetails,usedbins,UPthresh]=find_UPDOWN_cutoff(spikes,ledVals,window,powerRatio,r,bi,stimcond)

% window={[0 0.5]; ...
%         [0.25 0.75]; ...
%         [0.5 1]; ...
%         [0.75 1.25]; ...
%         [1 1.5]; ...
%         [1.25 1.75]; ...
%         [1.5 2]; ...
%         [1.75 2.25]; ...
%         [2 2.5]; ...
%         [2.25 2.75]; ...
%         [2.5 3]};
% led={[1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]};
% window={[1 1.5]+0.05; ...
%         [1.1 1.6]+0.05; ...
%         [1.2 1.7]+0.05; ...
%         [1.3 1.8]+0.05};
% led={[2 7]; ...
%         [2 7]; ...
%         [2 7]; ...
%         [2 7]};
% led={nan; ...
%         nan; ...
%         nan; ...
%         nan};

% stimcond=1:9;
spikes=filtspikes(spikes,0,'stimcond',stimcond);
windowWidth=0.250;
windowOffset=0.250;
% windowWidth=0.125;
% windowOffset=0.125;
% windowWidth=0.05;
% windowOffset=0.05;
% windowTimewindow=[1.25 1.325];
windowTimewindow=[0 5];
% ledVals=[5];
% ledVals=nan;
bins=65;

a=windowTimewindow(1);
b=windowTimewindow(2);
i=1;
while a<b
    if a+windowWidth>b
        break
    end
    window{i}=[a a+windowWidth];
    led{i}=ledVals;
    a=a+windowOffset;
    i=i+1;
end

alln=[];
for i=1:length(led)
    if isnan(ledVals)
        [m,s,n]=calcMeanAndStdDuringWindow(spikes,window{i});
    else
        [m,s,n]=calcMeanAndStdDuringWindow(filtspikes(spikes,0,'led',led{i}),window{i});
%         [m,s,n]=calcDifferenceAcrossWindows(filtspikes(spikes,0,'led',led{i}),window{i},[1.35 1.6]);
%         n=n(n>0);
%         new_n=[];
%         for j=1:100
%             ran=randsample(length(n),floor(length(n)*0.75));
%             new_n=[new_n; n(ran)];
%         end
%         n=new_n;
    end
%     ruse = 75 + 50.*randn(length(n),1);
    if ~isempty(r)
        % Make CDF from PDF
        % PDF is in r
        cdf=zeros(1,length(r.x));
        for j=1:length(r.x)
            if j==1
                cdf(j)=r.y(j);
            else
                cdf(j)=sum(r.y(1:j));
            end
        end
%         figure(); plot(cdf);
        % Draw samples from r using CDF
        b=max(cdf);
        a=0;
        samps=a+(b-a).*rand(length(n),1);
        rsamps=zeros(length(samps),1);
        for k=1:length(samps)
            if samps(k)==0
                rsamps(k)=0;
            else
                ind=find(cdf<samps(k),1,'last');
                rsamps(k)=r.x(ind);
            end
        end
%         figure(); 
%         [n_test,x_test]=hist(rsamps,50);
%         plot(x_test,n_test,'Color','k');
        % Add R to N
        alln=[alln; n+rsamps];
    else
        alln=[alln; n];
    end
    disp(i);
end

% bins=-5:10:800;
if ~isempty(bi)
    bins=bi;
end
[heights,centers]=hist(alln,bins);
usedbins=centers;
% [heights,centers]=histnorm(alln,5:5:500);
centers=[centers(1)-(centers(2)-centers(1)) centers centers(end)+(centers(2)-centers(1))];
heights=[0 heights 0];
figure(); 
hist(alln,bins);
% histnorm(alln,bins);
hold on;
plot(centers,heights);

% figure(); 
% histnorm(alln,bins);
% hold on;
% [fitparams,fitci]=fitPoissonPlusGaussian(alln);
% fitx=centers(1)-(centers(2)-centers(1)):1:centers(end)+(centers(2)-centers(1));
% y1=fitparams.a.*gampdf(fitx,fitparams.k,fitparams.theta);
% y2=fitparams.b.*normpdf(fitx,fitparams.mu,fitparams.sigma);
% y3=y1+y2;
% plot(fitx,y1,'Color','r');
% plot(fitx,y2,'Color','g');
% plot(fitx,y3,'Color','y');
% fitted=fitparams;

fitted=fit(centers',heights','gauss2');
% fitx=centers(1)-(centers(2)-centers(1)):0.1:centers(end)+(centers(2)-centers(1));
fitx=0:0.1:centers(end)+(centers(2)-centers(1));
y1=fitted.a1.*exp(-((fitx-fitted.b1)./fitted.c1).^2);
plot(fitx,y1,'Color','r');
y2=fitted.a2.*exp(-((fitx-fitted.b2)./fitted.c2).^2);
plot(fitx,y2,'Color','g');
y3=y1+y2;
plot(fitx,y3,'Color','y');
% interse=find(y2>y1,1);
interse=length(y2)-find(y2(end:-1:1)<y1(end:-1:1),1);
disp('x value at intersection');
disp(fitx(interse));
MUAthresh=fitx(interse);
histdetails.centers=centers;
histdetails.heights=heights;

UPthresh=3.694; % just initializing
stimDuration=5;
allpr=[];
if ~isempty(powerRatio)
    [a,uninds]=unique(spikes.trials);
    ledForSweeps=spikes.led(uninds);
    pr_xpoints=0:stimDuration/length(powerRatio{1}):stimDuration-(stimDuration/length(powerRatio{1}));
    for i=1:length(led)
        if isnan(ledVals)
            n=getAveragePR(powerRatio,window{i},pr_xpoints);
        else
            n=getAveragePR(powerRatio(ismember(ledForSweeps,led{i})),window{i},pr_xpoints);
        end
        allpr=[allpr; n];
    end

%     MUAthresh=30.3;
    muaUPs=alln>MUAthresh;
    threshSteps=0:0.1:10;
    n_errors=zeros(length(threshSteps),1);
    for i=1:length(threshSteps)
        UPthresh=threshSteps(i);
        prUPs=allpr>UPthresh;
        n_errors(i)=sum(muaUPs & ~prUPs) + sum(~muaUPs & prUPs);
    end
    [~,minerr]=min(n_errors);
    UPthresh=threshSteps(minerr);
    prUPs=allpr>UPthresh;
    disp('Best powerRatio UPthresh');
    disp(UPthresh);
    disp('Number of MUA UP states correctly classified as UP');
    disp(sum(muaUPs & prUPs));
    disp('Number of MUA UP states incorrectly classified as DOWN');
    disp(sum(muaUPs & ~prUPs));
    disp('Number of MUA DOWN states correctly classified as DOWN');
    disp(sum(~muaUPs & ~prUPs));
    disp('Number of MUA DOWN states incorrectly  classified as UP');
    disp(sum(~muaUPs & prUPs));
    disp('Total # errors');
    disp(sum(muaUPs & ~prUPs) + sum(~muaUPs & prUPs));
    
    figure(); 
    scatter(alln,allpr);
    hold on; 
    line([0 100],[UPthresh UPthresh]);
    scatter(alln(prUPs),allpr(prUPs),[],'k');
end

end    
function n=getAveragePR(powerRatio,window,pr_xpoints)
    
n=zeros(length(powerRatio),1);
for i=1:length(powerRatio)
    temp=powerRatio{i};
    n(i)=mean(temp(pr_xpoints>=window(1) & pr_xpoints<=window(2)));
end

end

function [fitparams,fitci]=fitPoissonPlusGaussian(data)

poiss_plus_gauss = @(x,a,b,mu,sigma,k,theta) (a.*gampdf(x,k,theta)+b.*normpdf(x,mu,sigma))./(a+b);

[n,x]=hist(data,50);

% start_mu=mean(data(data>mean(data)));
start_mu=75;
start_sigma=std(data(data>mean(data)));
start_a=1;
start_b=0.25;
start_k=2;
start_theta=mean(data)/2;
starts=[start_a start_b start_mu start_sigma start_k start_theta];
lowers=[0 0 0 0 0 0];

[phat,pci]=mle(data,'pdf',poiss_plus_gauss,'start',starts,'lowerbound',lowers);

fitparams.a=phat(1);
fitparams.b=phat(2);
fitparams.mu=phat(3);
fitparams.sigma=phat(4);
fitparams.k=phat(5);
fitparams.theta=phat(6);

fitci.a=pci(1);
fitci.b=pci(2);
fitci.mu=pci(3);
fitci.sigma=pci(4);
fitci.k=pci(5);
fitci.theta=pci(6);

end
         
