function measureMUASpontStimInteraction(spikes,fileInd)
% This function only works with MUA activity, not single units

spikes=filtspikes(spikes,0,'fileInd',fileInd);
useStimcond=1:8;
useLEDcond=[1 3 5 7];
useSpont_Stimcond=9;
stimOnset=1.05; % in s from start of trial
ledOnset=1.1;

% Get the normal progression of cortical activity from the moment of stim.
% onset during spontaneous activity
windowWidth=0.05; % not so much greater than 500 ms
if windowWidth>(ledOnset-stimOnset)
    disp('your predicted variable window overlaps with your independent window');
end
windowOffset=0.05;
spontStateWindow=0.05; % should be greater than windowWidth
% windowTimewindow=[1-spontStateWindow ledOnset];
windowTimewindow=[1.0 ledOnset];
ledVals=useLEDcond;

a=windowTimewindow(1);
b=windowTimewindow(2);
i=1;
while a<b
    if a+windowWidth>b
        break
    end
    window{i}=[a a+windowWidth];
    led{i}=ledVals;
    a=a+windowOffset;
    i=i+1;
end
stepbystep=cell(length(led),1);
av_stepbystep=cell(length(led),1);
for i=1:length(led)
    if isnan(ledVals)
        [m,s,n]=calcMeanAndStdDuringWindow(filtspikes(spikes,0,'stimcond',useSpont_Stimcond),window{i});
    else
        [m,s,n]=calcMeanAndStdDuringWindow(filtspikes(spikes,0,'led',led{i},'stimcond',useSpont_Stimcond),window{i});
    end
    stepbystep{i}=n;
    av_stepbystep{i}=mean(n);
    disp(i);
end
xpoints=zeros(length(window),1);
ypoints=zeros(length(window),1);
ypoints_trials=zeros(length(window),length(stepbystep{1}));
canuse_xpoints=zeros(length(window),1);
for i=1:length(window)
    xpoints(i)=mean(window{i});
    ypoints(i)=av_stepbystep{i};
    ypoints_trials(i,:)=stepbystep{i};
    temp=window{i};
    canuse_xpoints(i)=temp(2)<=stimOnset;
end
% ypoints_trials=ypoints_trials(:,(ypoints_trials(1,:)>200) );
disp([xpoints canuse_xpoints]);
figure(); 
scatter(ypoints_trials(1,:),ypoints_trials(end,:));
% figure();
% plot(xpoints,ypoints);
% figure(); 
% plot(xpoints,ypoints_trials);
% Multiple linear regression to predict spont. state at ledOnset from
% spont. state preceding stimulus onset
% include interaction terms?
% choose 50% of data for training model
set=1:size(ypoints_trials,2);
testset=set(1:4:end);
trainingset=set(~ismember(set,testset));
x=[];
for i=1:size(ypoints_trials(logical(canuse_xpoints),:),1)
    for j=i+1:size(ypoints_trials(logical(canuse_xpoints),:),1)
        x=[x (ypoints_trials(i,trainingset).*ypoints_trials(j,trainingset))'];
    end
end
% X=[ones(size(ypoints_trials,2),1) ypoints_trials(1:end-1,:)' x];
% X=[ones(size(ypoints_trials,2),1) ypoints_trials(1:end-1,:)'];
% X=[ones(size(ypoints_trials,2),1) ypoints_trials(logical(canuse_xpoints),:)'];
% X=[ones(size(ypoints_trials(:,1:2:end),2),1) ypoints_trials(logical(canuse_xpoints),1:2:end)' x];

X=[ones(size(ypoints_trials(:,trainingset),2),1) ypoints_trials(logical(canuse_xpoints),trainingset)'];
% X=[ones(size(ypoints_trials(:,trainingset),2),1) ypoints_trials(logical(canuse_xpoints),trainingset)' x];
[b,bint,r,rint,stats]=regress(ypoints_trials(end,trainingset)',X,0.05);
disp(b);
% figure(); 
% x1=min(X,1):0.001:max(X,1);
% x2=[];
% for i=2:size(X,2)
%     x2=[x2; min(X,i):0.001:max(X,2)];
% end
% yfit=b(1)*x1+b(2:end)*x2;
predictedy=zeros(length(ypoints_trials(end,:)),1);
actualy=ypoints_trials(end,:);
for i=1:size(ypoints_trials,2)
%     predictedy(i)=b(1)*1+sum(b(2:end).*ypoints_trials(logical(canuse_xpoints),i));
    currb=2+sum(canuse_xpoints);
    interactionTerms=0;
%     for j=1:size(ypoints_trials(logical(canuse_xpoints),:),1)
%         for k=j+1:size(ypoints_trials(logical(canuse_xpoints),:),1)
%             interactionTerms=interactionTerms+b(currb)*ypoints_trials(j,i)*ypoints_trials(k,i);
%             currb=currb+1;
%         end
%     end
%     predictedy(i)=b(1)*1+sum(b(2:1+sum(canuse_xpoints)).*ypoints_trials(logical(canuse_xpoints),i))+interactionTerms;
    predictedy(i)=b(1)*1+sum(b(2:1+sum(canuse_xpoints)).*ypoints_trials(logical(canuse_xpoints),i));
end
figure(); 
scatter(actualy(trainingset),predictedy(trainingset),[],'b');
hold on;
scatter(actualy(testset),predictedy(testset),[],'r');
disp('sum of test set residuals');
disp(sum(abs(actualy(testset)-predictedy(testset)'))/length(testset));
% miny=min(ypoints_trials(1,:));
% maxy=max(ypoints_trials(1,:));
% ybins=miny:(maxy-miny)/3:maxy;
% for i=1:length(ybins)-1
%    figure(); 
%    plot(xpoints,ypoints_trials(:,(ypoints_trials(1,:)>=ybins(i))&(ypoints_trials(1,:)<=ybins(i+1))));
% end
% figure(); 
% scatter(ypoints_trials(1,:),ypoints_trials(end,:));