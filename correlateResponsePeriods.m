function [correlations2,pvals,rlows,ruppers]=correlateResponsePeriods(data, meanOrSumOverPeriod, Fs, baselinePeriod, responsePeriods, useTheseRowsOfData, alignToBaseline)
% A figure will be made for each response period showing its correlation
% with the other response periods
% Also, all correlation data = a correlation matrix will be displayed as a 
% heat map, showing the correlation of each response period with all other
% response periods
% INPUT PARAMETERS
% data                  a matrix containing the response data
%                       each row should be a trial
%                       each column should be a different time point
%                       (sample)
% meanOrSumOverPeriod   a string indicating whether responses across trials
%                       (for each response period) should be summed or
%                       averaged
%                       'mean' -> average responses across trials
%                       'sum' -> sum across trials
% Fs                    sampling rate of data in Hz
% baselinePeriod        a two-element vector indicating the start time
%                       (relative to the beginning of a trial) for the 
%                       baseline period and the end time of the baseline
%                       period
%                       baselinePeriod(1) = time in seconds when baseline
%                       begins
%                       baselinePeriod(2) = time in seconds when baseline
%                       period ends
%                       * baselinePeriod may optionally be used to align
%                       all trials 
% responsePeriods       an n-by-2 array of times (in seconds relative to
%                       the start of a trial) indicating how to parse the
%                       responses into "periods", which are just consistent
%                       chunks of time across trials
%                       e.g., responsePeriods(1,1) = time in seconds when
%                       response period 1 begins
%                       responsePeriods(1,2) = time in seconds when
%                       response period 1 ends
%                       NOTE: the response periods and the baseline period
%                       are distinct (they are treated separately, as the
%                       baseline period may be used to align all trials)
% useTheseRowsOfData    a vector of the row indices of data to include in
%                       the analysis; all other rows of data will be
%                       ignored
% alignToBaseline       if alignToBaseline == 1, align all trials
%                       by subtracting for each trial its average during
%                       the baseline period
% OUTPUT
% correlations2         a correlation matrix giving the correlation of each
%                       response period with every other NORMALIZED by the 
%                       correlation of the response period with the baseline period 
%                       that is, correlations2 is n-by-n array where n is
%                       the number of response periods, and each element of
%                       correlations2 (i,j) is the correlation of the ith response period 
%                       with the jth response period MINUS the correlation of the ith 
%                       response period with the baseline period
% pvals                 a matrix corresponding to the correlation matrix,
%                       giving the pvals of the correlations
%                       where pvals(i,j) is the pval of the correlation of
%                       the ith response period with the jth response
%                       period
% rlows                 a matrix of the 95% confidence interval lower values 
%                       of the correlation coefficients 
% ruppers               a matrix of the 95% confidence interval upper
%                       values of the correlation coefficients

% Only use the specified trials
data=data(useTheseRowsOfData,:);

totalTrialLength=(1/Fs)*size(data,2);
baselineInds=floor(baselinePeriod(1)*Fs)+1:floor(baselinePeriod(2)*Fs);

% Align trials to baseline if alignToBaseline==1
if alignToBaseline
    data=data-mean(data(:,baselineInds),2)*ones(1,size(data,2));
end

% Get response periods in terms of columns of data
if size(responsePeriods,2)~=2
    disp('responsePeriods must be an n-by-2 matrix');
    return
end
responsePeriodsInInds=zeros(size(responsePeriods));
for i=1:size(responsePeriods,1)
    responsePeriodsInInds(i,:)=[floor(responsePeriods(i,1)*Fs)+1,floor(responsePeriods(i,2)*Fs)];
end
% Calculate responses for each period for each trial
if strcmp(meanOrSumOverPeriod,'mean')
    baselineResponses=mean(data(:,baselineInds),2);
else
    baselineResponses=sum(data(:,baselineInds),2);
end
responses=zeros(size(data,1),size(responsePeriods,1));
for i=1:size(responsePeriods,1)
    if strcmp(meanOrSumOverPeriod,'mean')
        responses(:,i)=mean(data(:,responsePeriodsInInds(i,1):responsePeriodsInInds(i,2)),2);
    else
        responses(:,i)=sum(data(:,responsePeriodsInInds(i,1):responsePeriodsInInds(i,2)),2);
    end 
end

figure;
scatter(responses(:,13),responses(:,20));

% Calculate correlations
correlations=zeros(size(responsePeriods,1)+1,size(responsePeriods,1)+1);
pvals=zeros(size(responsePeriods,1)+1,size(responsePeriods,1)+1);
rlows=zeros(size(responsePeriods,1)+1,size(responsePeriods,1)+1);
ruppers=zeros(size(responsePeriods,1)+1,size(responsePeriods,1)+1);
for i=0:size(responsePeriods,1)
    for j=i:size(responsePeriods,1)
        if i==j % response period's correlation with itself
            r=1;
            p=0;
            rlo=1;
            rup=1;
            correlations(i+1,j+1)=r;
            pvals(i+1,j+1)=p;
            rlows(i+1,j+1)=rlo;
            ruppers(i+1,j+1)=rup;
        else
            if i==0
                [r,p,rlo,rup]=corrcoef(baselineResponses,responses(:,j));
            else
                [r,p,rlo,rup]=corrcoef(responses(:,i),responses(:,j));
            end
            correlations(i+1,j+1)=r(1,2);
            pvals(i+1,j+1)=p(1,2);
            rlows(i+1,j+1)=rlo(1,2);
            ruppers(i+1,j+1)=rup(1,2);
            % Make matrices symmetric
            correlations(j+1,i+1)=r(1,2);
            pvals(j+1,i+1)=p(1,2);
            rlows(j+1,i+1)=rlo(1,2);
            ruppers(j+1,i+1)=rup(1,2);
        end
    end
end

% Plot each response period's correlation with baseline period
% figure(1);
% title('Correlation of each response period with baseline');
% bar(1:size(correlations,1),correlations(1,:));    
% 
% % Plot correlation results
% for i=2:size(correlations,1)
%     figure(i);
%     title(['Correlations of Response Period ' num2str(i-1) ' with Other Periods']);
%     bar(1:size(correlations,1),correlations(i,:));
% end

% Plot correlation results with respect to response period correlation with
% baseline
correlations2=zeros(size(correlations));
% for i=2:size(correlations,1)
%     correlations2(i,:)=correlations(i,:)-correlations(1,:);
% end
for i=2:size(correlations,1)
    correlations2(i,:)=correlations(i,:)-correlations(i,1);
end
% for i=2:size(correlations2,1)
%     figure(i+size(correlations2,1));
%     title(['Corr. of Response Period ' num2str(i-1) ' with Other Periods, Normalized to Its Corr. with Baseline']);
%     %bar(1:size(correlations2,1),correlations2(i,:));
%     subplot(3,1,1);
%     plot([mean(baselinePeriod), mean(responsePeriods,2)'],correlations2(i,:));
%     subplot(3,1,2);
%     hmap_vector=zeros(1,size(data,2));
%     hmap_vector(baselineInds)=0;
%     for j=1:size(responsePeriodsInInds,1)
%         hmap_vector(responsePeriodsInInds(j,1):responsePeriodsInInds(j,2))=correlations2(i,j+1);
%     end
%     heatmap(hmap_vector,[],[],'','TextColor','xor','FontSize',16,'Colorbar',{'EastOutside'});
%     colormap jet;
%     subplot(3,1,3);
%     plot(0:totalTrialLength/(size(data,2)-1):totalTrialLength,mean(data,1),'Color','k');
%     hold on;
%     startTime=responsePeriods(i-1,1);
%     endTime=responsePeriods(i-1,2);
%     startInds=responsePeriodsInInds(i-1,1);
%     endInds=responsePeriodsInInds(i-1,2);
%     plot(startTime:(endTime-startTime)/(length(startInds:endInds)-1):endTime,mean(data(:,startInds:endInds),1),'Color','r','Linewidth',1.5);
% %     for j=0:size(responsePeriods)
% %         if j==0
% %             currInds=baselineInds;
% %         else
% %             currInds=responsePeriodsInInds(j,1):responsePeriodInInds(j,2);
% %         end
% %         c=
% %         plot(responsePeriods(j,1):responsePeriods(j,2)/length(currInds)-1:responsePeriods(j,2),mean(data(:,currInds),1),'Color','');
% end

% Make colormap of correlations
% figure();
% heatmap(correlations,[],[],'','TextColor','xor','FontSize',16,'Colorbar',{'EastOutside'});
% colormap jet;
% Heatmap of p-values
figure();
pvals(pvals>0.05)=0.1;
heatmap(pvals,[],[],'','TextColor','xor','FontSize',16,'Colorbar',{'EastOutside'});
colormap jet;
% Heatmap of significance
% pvals(pvals>0.05)=0;
% pvals(pvals<=0.05)=1;
% figure();
% heatmap(pvals,[],[],'','TextColor','xor','FontSize',16,'Colorbar',{'EastOutside'});
% colormap jet;

% Make colormap of correlations but scale to response periods and show mean
% signal
% hmap_matrix=zeros(size(data,2),size(data,2));
% for i=0:size(responsePeriodsInInds,1)
%     for j=i:size(responsePeriodsInInds,1)
%         if i==0 && j==0
%             hmap_matrix(baselineInds,baselineInds)=correlations(1,1);
%         elseif i==0
%             hmap_matrix(baselineInds,responsePeriodsInInds(j,1):responsePeriodsInInds(j,2))=correlations(1,j+1);
%             hmap_matrix(responsePeriodsInInds(j,1):responsePeriodsInInds(j,2),baselineInds)=correlations(1,j+1);
%         else
%             hmap_matrix(responsePeriodsInInds(i,1):responsePeriodsInInds(i,2),responsePeriodsInInds(j,1):responsePeriodsInInds(j,2))=correlations(i+1,j+1);
%             hmap_matrix(responsePeriodsInInds(j,1):responsePeriodsInInds(j,2),responsePeriodsInInds(i,1):responsePeriodsInInds(i,2))=correlations(i+1,j+1);
%         end
%     end
% end
hmap_matrix=zeros(size(data,2),size(data,2));
for i=1:size(responsePeriodsInInds,1)
    for j=i:size(responsePeriodsInInds,1)
        hmap_matrix(responsePeriodsInInds(i,1):responsePeriodsInInds(i,2),responsePeriodsInInds(j,1):responsePeriodsInInds(j,2))=correlations(i+1,j+1);
        hmap_matrix(responsePeriodsInInds(j,1):responsePeriodsInInds(j,2),responsePeriodsInInds(i,1):responsePeriodsInInds(i,2))=correlations(i+1,j+1);
    end
end
figure();
subplot(2,1,1);
%heatmap(hmap_matrix,[],[],'','TextColor','xor','FontSize',16,'Colorbar',{'EastOutside'});
heatmap(hmap_matrix,[],[],'','TextColor','xor','FontSize',16);
colormap jet;
subplot(2,1,2);
plot(0:totalTrialLength/(size(data,2)-1):totalTrialLength,mean(data,1),'Color','k');
axis tight;

disp('Correlation of 0.72-0.9 s period with 1.3-1.5 s period');
disp(correlations(8,11));
disp(pvals(8,11));
disp('Correlation of baseline period with 1.3-1.5 s period');
disp(correlations(1,11));
disp(pvals(1,11));
disp('Correlation of 0.575-0.59 s period with 1.3-1.5 s period');
disp(correlations(4,11));
disp(pvals(4,11));