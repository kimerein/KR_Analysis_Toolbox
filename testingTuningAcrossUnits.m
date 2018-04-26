function [control_OSIs LED_OSIs]=testingTuningAcrossUnits(tunedFRs)

% Calculate OSIs for control and LED on
control_OSIs=[];
LED_OSIs=[];
for i=1:length(tunedFRs)
    t=tunedFRs{i};
    if isempty(t)
        continue
    end
    [pref_fr,pref_orient]=max(t(:,1));
    [nonpref_fr,nonpref_orient]=min(t(:,1));
    control_OSIs=[control_OSIs; pref_fr/(pref_fr+nonpref_fr)];
    [pref_fr,pref_orient]=max(t(:,2));
    [nonpref_fr,nonpref_orient]=min(t(:,2));
    LED_OSIs=[LED_OSIs; pref_fr/(pref_fr+nonpref_fr)];
end
figure(); 
plot([control_OSIs LED_OSIs]');


numEmpty=0;
for i=1:length(tunedFRs)
    % Find pref. orient.
    t=tunedFRs{i};
    if isempty(t)
        numEmpty=numEmpty+1;
        continue
    end
    [pref_fr,pref_orient]=max(t(:,1));
    shift_t=zeros(size(t));
    k=1;
    for j=pref_orient:length(t(:,1))
        shift_t(k,:)=t(j,:);
        k=k+1;
    end
    for j=1:pref_orient-1
        shift_t(k,:)=t(j,:);
        k=k+1;
    end
    T{i}=shift_t;
end

mean_tunedFRs=zeros(size(T{1}));
for i=1:length(T)
    if isempty(T{i})
        continue
    end
    mean_tunedFRs=mean_tunedFRs+T{i};
end
mean_tunedFRs=mean_tunedFRs./(length(tunedFRs)-numEmpty);

figure(); 
% plot([1:8]',[mean_tunedFRs(5:8,1); mean_tunedFRs(1:4,1)],'Color','k');
% hold on;
% plot([1:8]',[mean_tunedFRs(5:8,2); mean_tunedFRs(1:4,2)],'Color','r');
plot([1:9]',[mean_tunedFRs(:,1); mean_tunedFRs(1,1)],'Color','k');
hold on;
plot([1:9]',[mean_tunedFRs(:,2); mean_tunedFRs(1,2)],'Color','r');

% for i=1:length(T)
%     if isempty(T{i})
%         continue
%     end
%     t=T{i};
%     [min_FR,min_ind]=min(t(:,1));
%     max_FR=t(1,1);
%     max_ind=1;
%         