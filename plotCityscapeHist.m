function plotCityscapeHist(data1,data2,bins)

normByIntegral=false;

if isscalar(bins)
    mi=nanmin([data1; data2]);
    ma=nanmax([data1; data2]);
    binwidth=(ma-mi)/bins;
    % center bins around 0
    binsBelowZero=0-(binwidth/2):-binwidth:mi;
    binsAboveZero=0+(binwidth/2):binwidth:ma;
    bins=[fliplr(binsBelowZero) binsAboveZero];
end

[n,x]=histcounts(data1,bins);
[n,x]=cityscape_hist(n,x);
figure();
if normByIntegral==true
    plot(x,n./nansum(n),'Color','b');
else
    plot(x,n,'Color','b');
end
if isempty(data2)
    return
end
hold on;
[n,x]=histcounts(data2,bins);
[n,x]=cityscape_hist(n,x);
if normByIntegral==true
    plot(x,n./nansum(n),'Color','r');
else
    plot(x,n,'Color','r');
end
legend({'data 1','data 2'});

end