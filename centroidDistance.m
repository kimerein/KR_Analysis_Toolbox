function [cents,withinGroupDistance,acrossGroupDistance]=centroidDistance(data)

% First column of data should be group classification
% Other columns are different dimensions of the data

groups=unique(data(:,1));
cents=nan(length(groups),size(data,2));
withinGroupDistance=nan(length(groups),1);
for i=1:length(groups)
    currGroup=groups(i);
    subData=data(data(:,1)==currGroup,:);
    cents(i,:)=sum(subData,1)/size(subData,1);
    sumDistance=0;
    for j=1:size(subData,1)
        di=subData(j,2:end)-cents(i,2:end);
        disq=di.^2;
        sumDistance=sumDistance+sqrt(sum(disq));
    end        
    withinGroupDistance(i)=sumDistance;
end

acrossGroupDistance=nan(length(groups),1);
for i=1:length(groups)
    ccent=cents(i,:);
    othercent=cents([1:i-1 i+1:end],:);
    sumDistance=0;
    for j=1:size(othercent,1)
        di=othercent(j,2:end)-ccent(2:end);
        disq=di.^2;
        sumDistance=sumDistance+sqrt(sum(disq));
    end
    acrossGroupDistance(i)=sumDistance;
end
