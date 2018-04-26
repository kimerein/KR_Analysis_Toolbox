function [pvals,final_pval]=compareFRHistograms_figure(spikes,useTheseAssigns,stimCond_black,ledCond_black,stimCond_red,ledCond_red,timeWindow_black,timeWindow_red,trials,nHistBins)

figure(); 
numSubplots=length(stimCond_black);
blackAvs=zeros(numSubplots,1);
redAvs=zeros(numSubplots,1);
pvals=zeros(numSubplots,1);
b_n=[];
r_n=[];
evSub_n=[];
for i=1:numSubplots
    [pvals(i),blackAvs(i),redAvs(i),blackPlot,redPlot,black_n,red_n]=compareFRHistograms_returnAxes(spikes,useTheseAssigns,stimCond_black{i},ledCond_black{i},stimCond_red{i},ledCond_red{i},timeWindow_black,timeWindow_red,trials,nHistBins);
    b_n=[b_n; black_n];
%     b_n=black_n;
    r_n=[r_n; red_n];
%     evSub_n=[evSub_n; red_n-black_n];
    h=subplot(numSubplots,1,i);
    plot(blackPlot.x,blackPlot.y,'Color','k');
    hold on; 
    plot(redPlot.x,redPlot.y,'Color','r');
%     [sub_n,sub_x]=hist(red_n-black_n,nHistBins);
%     plot(sub_x,sub_n,'Color','blue');
%     set(h,'XLim',[-400 400]);
    if i==numSubplots
        xlabel('Firing Rate (Hz)');
        ylabel('Count');
    end
end

figure(); 
[all_n,all_x]=histnorm(b_n,nHistBins);
plot(all_x,all_n,'Color','black');
[all_n,all_x]=histnorm(r_n,nHistBins);
hold on; 
plot(all_x,all_n,'Color','red');
% [all_n,all_x]=histnorm(evSub_n,nHistBins);
% plot(all_x,all_n,'Color','blue');

[final_pval]=ranksum(b_n,r_n);