function [control_OSIs LED_OSIs]=circularVariance(tunedFRs)

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