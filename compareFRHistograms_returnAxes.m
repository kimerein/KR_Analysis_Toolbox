function [p,black_m,red_m,blackPlot,redPlot,black_n,red_n]=compareFRHistograms_returnAxes(spikes,useTheseAssigns,stimCond_black,ledCond_black,stimCond_red,ledCond_red,timeWindow_black,timeWindow_red,trials,nHistBins)

if isempty(useTheseAssigns)
else
    spikes=filtspikes(spikes,0,'assigns',useTheseAssigns);
end

if isempty(trials)
else
    spikes=filtspikes(spikes,0,'trials',trials);
end

% Get FR distribution for black conditions
black_spikes=filtspikes(spikes,0,'stimcond',stimCond_black,'led',ledCond_black);
[black_m,black_s,black_n]=calcMeanAndStdDuringWindow(black_spikes,timeWindow_black);

% Get FR distribution for red conditions
red_spikes=filtspikes(spikes,0,'stimcond',stimCond_red,'led',ledCond_red);
[red_m,red_s,red_n]=calcMeanAndStdDuringWindow(red_spikes,timeWindow_red);

% Plot figure
% figure(); 
[b_n,b_x]=hist(black_n,nHistBins);
[r_n,r_x]=hist(red_n,nHistBins);
% plot(b_x,b_n,'Color','black');
% hold on; 
% plot(r_x,r_n,'Color','red');
blackPlot.x=b_x;
blackPlot.y=b_n;
redPlot.x=r_x;
redPlot.y=r_n;

[p,h,stats]=ranksum(black_n,red_n);