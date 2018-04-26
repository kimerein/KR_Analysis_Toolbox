function [h_diffs,p_diffs]=putLayerDataTogetherAcrossAnimals(layerData,dataBlock)

n_noLED=zeros(length(layerData),length(layerData{1}));
n_LED=zeros(length(layerData),length(layerData{1}));

for i=1:length(layerData)
    curr=layerData{i};
    for j=1:length(layerData{i})
        n_noLED(i,j)=mean(curr{j}.FRs_noLED{dataBlock});
        n_LED(i,j)=mean(curr{j}.FRs_LED{dataBlock});
%         n_noLED(i,j)=mean(curr.FRs_noLED{dataBlock});
%         n_LED(i,j)=mean(curr.FRs_LED{dataBlock});
    end
end

figure(); 
m_noLED=mean(n_noLED,1);
m_LED=mean(n_LED,1);
std_noLED=std(n_noLED,[],1);
std_LED=std(n_LED,[],1);
errorbar(1:length(layerData{1}),m_noLED,std_noLED,'Color','k');
hold on;
errorbar(1:length(layerData{1}),m_LED,std_LED,'Color','r');

backup_n_noLED=n_noLED;
backup_n_LED=n_LED;

for i=1:size(n_noLED,1)
    minVal=min(n_noLED(i,:));
%     minVal=min(n_LED(i,:));
    n_noLED(i,:)=n_noLED(i,:)-minVal;
    n_LED(i,:)=n_LED(i,:)-minVal;
    maxVal=max(n_noLED(i,:));
%     maxVal=max(n_LED(i,:));
    n_noLED(i,:)=n_noLED(i,:)*(1/maxVal);
    n_LED(i,:)=n_LED(i,:)*(1/maxVal);
end
figure(); 
m_noLED=mean(n_noLED,1);
m_LED=mean(n_LED,1);
std_noLED=std(n_noLED,[],1);
std_LED=std(n_LED,[],1);
errorbar(1:length(layerData{1}),m_noLED,std_noLED,'Color','k');
hold on;
errorbar(1:length(layerData{1}),m_LED,std_LED,'Color','r');

% n_noLED=backup_n_noLED;
% n_LED=backup_n_LED;

% figure(); 
% plot(1:length(layerData{1}),n_noLED(5,:),'Color','k');
% hold on;
% plot(1:length(layerData{1}),n_LED(5,:),'Color','r');

diffs=zeros(size(n_noLED));
for i=1:size(n_noLED,1)
    diffs(i,:)=n_noLED(i,:)-n_LED(i,:);
end
figure(); 
m_diffs=mean(diffs,1);
std_diffs=std(diffs,[],1);
errorbar(1:length(layerData{1}),m_diffs,std_diffs);

h_diffs=zeros(1,size(n_noLED,2));
p_diffs=zeros(1,size(n_noLED,2));
for i=1:size(n_noLED,2)
    [h_diffs(i),p_diffs(i)]=ttest(n_noLED(:,i),n_LED(:,i),0.05,'both');
end
% [h_matrix,p_matrix]=makeSignificanceMatrix(diffs,0.05,'both');
% figure(); 
% heatmap(h_matrix);