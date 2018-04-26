function [diags1,diags2,best_pvals,diags1scale,diags2scale]=alignVolumes(diags1,diags2,baseVol1,baseVol2,bestS)

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];


[diags1,vol1.mid,vol1.lower,vol1.upper]=plotFreqResponse_upperLimit(diags1,[],[],50,30,70,0,[]);
[diags2,vol2.mid,vol2.lower,vol2.upper]=plotFreqResponse_upperLimit(diags2,[],[],50,30,70,0,[]);

vol1.upper=vol1.upper-baseVol1;
vol1.mid=vol1.mid-baseVol1;
vol1.lower=vol1.lower-baseVol1;
vol2.upper=vol2.upper-baseVol2;
vol2.mid=vol2.mid-baseVol2;
vol2.lower=vol2.lower-baseVol2;
diags1=diags1-baseVol1;
diags2=diags2-baseVol2;

temp1=vol1.mid;
temp2=vol2.mid;
diags1=diags1./mean(temp1);
diags1scale=1/mean(temp1);
diags2=diags2./mean(temp2);
diags2scale=1/mean(temp2);
vol1.upper=vol1.upper./mean(temp1);
vol1.mid=vol1.mid./mean(temp1);
vol1.lower=vol1.lower./mean(temp1);
vol2.upper=vol2.upper./mean(temp2);
vol2.mid=vol2.mid./mean(temp2);
vol2.lower=vol2.lower./mean(temp2);

% minscale=min(vol1.mid(1:10))./mean(temp2);
% maxscale=max(vol1.upper(1:10))./mean(temp2);
% tryscales=linspace(minscale,maxscale,100);
tryscales=0.05:0.1:10;

% Get summed p-val
summed_pvals=nan(1,length(tryscales));
summed_binary_pvals=nan(1,length(tryscales));
for i=1:length(tryscales)
    curr_diags1=diags1.*tryscales(i);
    running_psum=0;
    running_binary_psum=0;
    for j=1:15
        currp=ranksum(curr_diags1(:,j),diags2(:,j));
        running_psum=running_psum+currp;
        running_binary_psum=running_binary_psum+(currp>0.05);
    end
    summed_pvals(i)=running_psum;
    summed_binary_pvals(i)=running_binary_psum;
end

% Maximize summed p-val
m=max(summed_binary_pvals);
useIndStep1=find(summed_binary_pvals==m);
[~,m2_ind]=max(summed_pvals(useIndStep1));
final_ind=useIndStep1(m2_ind);
bestScale=tryscales(final_ind);
if ~isempty(bestS)
    bestScale=bestS;
end
diags1=diags1.*bestScale;
diags1scale=diags1scale*bestScale;
    
best_pvals=nan(1,15);
for j=1:15
    currp=ranksum(diags1(:,j),diags2(:,j));
    best_pvals(j)=currp;
end
disp(best_pvals);
