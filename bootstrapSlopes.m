function [meanSlope,meanIntercept,SEMslope,SEMintercept,avSlopes,avIncps,bootIndExptFits]=bootstrapSlopes(fEv,fSpont)

% fEv = fractional evoked suppression
% in format
% x rows by y columns
% where columns are different values for different LED conds
% and rows are different experiments or trials
% 
% fSpont = fractional spontaneous suppression
% must match size of fEv
% and map onto fEv point-by-point

% Do each experiment or row individually
% Bootstrap slope values for that experiment
% Get average slope for that experiment
% And deviation

n=100; % number of bootstrap trials per experiment
g=fittype('m*x','coefficients','m','independent','x');
for j=1:n
    % Get random permutation matrix with replacement
    for i=1:size(fEv,1)
        pEv(i,:)=randsample(size(fEv,2),size(fEv,2),true);
        while all(pEv(i,:)==pEv(i,1))
            pEv(i,:)=randsample(size(fEv,2),size(fEv,2),true);
        end
    end
    for i=1:size(fEv,1)
        thesefEvs(i,:)=fEv(i,pEv(i,:));
        thesefSponts(i,:)=fSpont(i,pEv(i,:));
    end
    % Calculate best-fit slope for each sample for each experiment
    for i=1:size(thesefEvs,1)
        % To have fit choose y-intercept
%         p=polyfit(thesefEvs(i,:),thesefSponts(i,:),1);
%         bestFits(i,:)=p;
%         yfit=polyval(p,thesefEvs(i,:));

        % To require a y-intercept=0
        fitobject=fit(thesefEvs(i,:)',thesefSponts(i,:)',g);
        bestFits(i,:)=[fitobject.m 0];
        yfit=fitobject.m.*thesefEvs(i,:);
      
        yresid=thesefSponts(i,:)-yfit;
        SSresid=sum(yresid.^2);
        SStotal=(length(thesefSponts(i,:))-1)*var(thesefSponts(i,:));
        rsq=1-SSresid/SStotal;
        rsqs(i)=rsq;
    end
    for i=1:size(thesefEvs,1)
        bootIndExptFits(j,i,:)=bestFits(i,:);
        bootIndExptRsqs(j,i)=rsqs(i);
    end
end
% Get average slope for each experiment
for i=1:size(fEv,1)
    avSlopes(i)=mean(bootIndExptFits(:,i,1));
    avIncps(i)=mean(bootIndExptFits(:,i,2));
end
% Average slopes across experiments and get S.E.M. slope
meanSlope=mean(avSlopes);
meanIntercept=mean(avIncps);
SEMslope=std(avSlopes);
SEMintercept=std(avIncps);