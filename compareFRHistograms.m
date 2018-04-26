function [p,black_m,red_m,black_s,red_s]=compareFRHistograms(spikes,useTheseAssigns,stimCond_black,ledCond_black,stimCond_red,ledCond_red,timeWindow_black,timeWindow_red,trials,nHistBins)

enforceSameBins=1;
normalizeHistograms=1;
useTheseBinCenters=[0 10 20 40 80 160 320 640 1280];
if ~isempty(useTheseBinCenters)
    nHistBins=useTheseBinCenters;
end

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
if normalizeHistograms==0
    figure(); 
    [b_n,b_x]=hist(black_n,nHistBins);
    if enforceSameBins==1
        [r_n,r_x]=hist(red_n,nHistBins);
    else
        [r_n,r_x]=hist(red_n,b_x);
    end
    semilogx(b_x,b_n,'Color','black');
    hold on;
    semilogx(r_x,r_n,'Color','red');
else
    figure(); 
    [b_n,b_x]=histnorm(black_n,nHistBins);
    if enforceSameBins==1
        [r_n,r_x]=histnorm(red_n,nHistBins);
    else
        [r_n,r_x]=histnorm(red_n,b_x);
    end
    semilogx(b_x,b_n,'Color','black');
    hold on;
    semilogx(r_x,r_n,'Color','red');
end

[p,h,stats]=ranksum(black_n,red_n);