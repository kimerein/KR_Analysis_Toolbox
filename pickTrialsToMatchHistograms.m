function [r]=pickTrialsToMatchHistograms(spikes,stimCond_black,ledCond_black,stimCond_red,ledCond_red,timeWindow,nHistBins)

% Get FR distribution for black conditions
black_spikes=filtspikes(spikes,0,'stimcond',stimCond_black,'led',ledCond_black);
[black_m,black_s,black_n]=calcMeanAndStdDuringWindow(black_spikes,timeWindow);

% Get FR distribution for red conditions
red_spikes=filtspikes(spikes,0,'stimcond',stimCond_red,'led',ledCond_red);
[red_m,red_s,red_n]=calcMeanAndStdDuringWindow(red_spikes,timeWindow);
redTrialInds=unique(red_spikes.trials);

% Plot figure before trial choice
figure(); 
[b_n,b_x]=hist(black_n,nHistBins);
[r_n,r_x]=hist(red_n,b_x);
plot(b_x,b_n,'Color','black');
hold on; 
plot(r_x,r_n,'Color','red');
[p,h,stats]=ranksum(black_n,red_n);

% Get trial by trial bins in histogram
histBinSize=b_x(2)-b_x(1);
binStarts=b_x(1)-histBinSize/2:histBinSize:b_x(end)-histBinSize/2;
binEnds=b_x(1)+histBinSize/2:histBinSize:b_x(end)+histBinSize/2;
blackTrialBins=zeros(1,length(black_n));
for i=1:length(black_n)
    for j=1:length(binStarts)
        if black_n(i)>=binStarts(j) & black_n(i)<binEnds(j)
            blackTrialBins(i)=j;
            break
        end
    end
end
redTrialBins=zeros(1,length(red_n));
for i=1:length(red_n)
    for j=1:length(binStarts)
        if red_n(i)>=binStarts(j) & red_n(i)<binEnds(j)
            redTrialBins(i)=j;
            break
        end
    end
end
    
rat=b_n./r_n;
% rat(isinf(rat))=b_n(isinf(rat))/(1/(sum(r_n)+1));
min_rat=min(rat);
rat(isinf(rat))=min_rat;
[max_rat,max_ind]=max(rat);
useTheseRedTrials=[];
for i=1:length(rat)
    if i==max_ind
        continue
    end
    nBlackTrials=b_n(i);
    useNRedTrials=ceil(nBlackTrials/max_rat);
    couldUseRedInds=find(redTrialBins==i);
    ncoulduse=length(couldUseRedInds);
    if ncoulduse==0
        continue
    end
    useTheseRedTrials=[useTheseRedTrials couldUseRedInds(randi(ncoulduse,[1 useNRedTrials]))];
end

% Plot figure after trial choice
figure(); 
[b_n,b_x]=hist(black_n,nHistBins);
[r_n,r_x]=hist(red_n(useTheseRedTrials),b_x);
plot(b_x,b_n,'Color','black');
hold on; 
plot(r_x,r_n,'Color','red');
plot(r_x,r_n*max_rat,'Color','c');
[p,h,stats]=ranksum(black_n,red_n(useTheseRedTrials));
disp(p);
r=redTrialInds(useTheseRedTrials);
[red2_m,red2_s,red2_n]=calcMeanAndStdDuringWindow(filtspikes(red_spikes,0,'trials',r),timeWindow);
figure(); 
[r_n,r_x]=hist(red_n(useTheseRedTrials),nHistBins);
[r2_n,r2_x]=hist(red2_n,r_x);
plot(r_x,r_n,'Color','r');
hold on; 
plot(r2_x,r2_n,'Color','green');