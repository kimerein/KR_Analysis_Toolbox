function [x,psth1,psth2,psth3]=doPredictionAnalysis(overall_l1,overall_l2,overall_l3,alphaYes)

[t1_1,x]=plotUnitResponseTransition(filtspikes(overall_l1,0,'alpha',alphaYes),unique(overall_l1.assigns),1,[1 2 3]);
[t2_1,x]=plotUnitResponseTransition(filtspikes(overall_l2,0,'alpha',alphaYes),unique(overall_l2.assigns),1,[1 2 3]);
[t3_1,x]=plotUnitResponseTransition(filtspikes(overall_l3,0,'alpha',alphaYes),unique(overall_l3.assigns),1,[1 2 3]);
close all;

[t1_2,x]=plotUnitResponseTransition(filtspikes(overall_l1,0,'alpha',alphaYes),unique(overall_l1.assigns),1,[4 5 6]);
[t2_2,x]=plotUnitResponseTransition(filtspikes(overall_l2,0,'alpha',alphaYes),unique(overall_l2.assigns),1,[4 5 6]);
[t3_2,x]=plotUnitResponseTransition(filtspikes(overall_l3,0,'alpha',alphaYes),unique(overall_l3.assigns),1,[4 5 6]);
close all;

[t1_3,x]=plotUnitResponseTransition(filtspikes(overall_l1,0,'alpha',alphaYes),unique(overall_l1.assigns),1,[7 8 9]);
[t2_3,x]=plotUnitResponseTransition(filtspikes(overall_l2,0,'alpha',alphaYes),unique(overall_l2.assigns),1,[7 8 9]);
[t3_3,x]=plotUnitResponseTransition(filtspikes(overall_l3,0,'alpha',alphaYes),unique(overall_l3.assigns),1,[7 8 9]);
close all;

[t1_4,x]=plotUnitResponseTransition(filtspikes(overall_l1,0,'alpha',alphaYes),unique(overall_l1.assigns),1,[10 11 12]);
[t2_4,x]=plotUnitResponseTransition(filtspikes(overall_l2,0,'alpha',alphaYes),unique(overall_l2.assigns),1,[10 11 12]);
[t3_4,x]=plotUnitResponseTransition(filtspikes(overall_l3,0,'alpha',alphaYes),unique(overall_l3.assigns),1,[10 11 12]);
close all;

[x,psth1,psth2,psth3]=plotPvaluePSTH(x,[t1_1; t1_2; t1_3; t1_4],[t2_1; t2_2; t2_3; t2_4],[t3_1; t3_2; t3_3; t3_4],50);