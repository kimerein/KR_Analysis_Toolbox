function [bestC,bestDiff,ranks1,ranks2]=alignCSDs(dataInt1,dataInt2,window1,window2)
% dataInt1 and dataInt2 must have the same number of rows

Fs=32000;
downSamped=10;
Fs=Fs/downSamped;
times=0:1/Fs:size(dataInt1,2)*(1/Fs);

for i=1:size(dataInt1,1)
    subset=dataInt1(i,times>=window1(1) & times<=window1(2));
    [sink(i),sinkPeakInd(i)]=findFirstLocalMin(subset);
end

[sorted,inds1]=sort(sinkPeakInd);
a=1:length(sorted);
ranks1=a(inds1);
ranks1=smooth(ranks1,5);
figure(); 
plot(ranks1,1:length(ranks1));

for i=1:size(dataInt2,1)
    subset=dataInt2(i,times>=window2(1) & times<=window2(2));
    [sink(i),sinkPeakInd(i)]=findFirstLocalMin(subset);
end

[sorted,inds2]=sort(sinkPeakInd);
a=1:length(sorted);
ranks2=a(inds2);
ranks2=smooth(ranks2,5);
figure(); 
plot(ranks2,1:length(ranks2));

[a,r1_ind]=min(ranks1);
[b,r2_ind]=min(ranks2);
bestC=r1_ind-r2_ind;
bestDiff=0;

% bestC=0;
% bestDiff=nan;
% % for c=0:length(ranks1)-1
% for c=0:80
%     runningDiff=0;
%     chsAligned=0;
%     for j=1:length(ranks1)
%         if j+c>length(ranks2)
%         else
%             chsAligned=chsAligned+1;
%             runningDiff=runningDiff+abs(ranks1(j)-ranks2(j+c));
%         end
%     end
%     runningDiff=runningDiff/chsAligned;
%     if isnan(bestDiff)
%         bestDiff=runningDiff;
%         bestC=c;
%     elseif runningDiff<bestDiff
%         bestDiff=runningDiff;
%         bestC=c;
%     end
% end

figure(); 
plot(ranks1,1:length(ranks1),'Color','k');
hold on; 
plot(ranks2,bestC+1:bestC+length(ranks2),'Color','b');
end

function [curr,lowestInd]=findFirstLocalMin(vec)

curr=nan;
lowestInd=nan;
for i=1:length(vec)
    if isnan(curr)
        curr=vec(i);
        lowestInd=i;
    elseif i==length(vec)-1
        curr=vec(i);
        lowestInd=i;
        break
    elseif vec(i)<vec(i-1) && vec(i)<vec(i+1)
        curr=vec(i);
        lowestInd=i;
        break
    end
end
end