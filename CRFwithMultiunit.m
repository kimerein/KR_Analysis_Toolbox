% No LED trials to show contrast response function to drifting grating

% newSpikes = makeTempField(spikes,'led',0.01);
% T1_noLEDspikes = filtspikes(newSpikes,0,'temp',1);

% For Stim. Block 17 used in 3/30/12, collapse across orientation
% Just separate by DG contrast
% Orientation was first variable

% T1
T1_con1=filtspikes(T1_noLEDspikes,0,'stimcond',[1 7 13]);
T1_con2=filtspikes(T1_noLEDspikes,0,'stimcond',[2 8 14]);
T1_con3=filtspikes(T1_noLEDspikes,0,'stimcond',[3 9 15]);
T1_con4=filtspikes(T1_noLEDspikes,0,'stimcond',[4 10 16]);
T1_con5=filtspikes(T1_noLEDspikes,0,'stimcond',[5 11 17]);
T1_con6=filtspikes(T1_noLEDspikes,0,'stimcond',[6 12 18]);

% T3
T2_con1=filtspikes(T2_noLEDspikes,0,'stimcond',[1 7 13]);
T2_con2=filtspikes(T2_noLEDspikes,0,'stimcond',[2 8 14]);
T2_con3=filtspikes(T2_noLEDspikes,0,'stimcond',[3 9 15]);
T2_con4=filtspikes(T2_noLEDspikes,0,'stimcond',[4 10 16]);
T2_con5=filtspikes(T2_noLEDspikes,0,'stimcond',[5 11 17]);
T2_con6=filtspikes(T2_noLEDspikes,0,'stimcond',[6 12 18]);

% T3
T3_con1=filtspikes(T3_noLEDspikes,0,'stimcond',[1 7 13]);
T3_con2=filtspikes(T3_noLEDspikes,0,'stimcond',[2 8 14]);
T3_con3=filtspikes(T3_noLEDspikes,0,'stimcond',[3 9 15]);
T3_con4=filtspikes(T3_noLEDspikes,0,'stimcond',[4 10 16]);
T3_con5=filtspikes(T3_noLEDspikes,0,'stimcond',[5 11 17]);
T3_con6=filtspikes(T3_noLEDspikes,0,'stimcond',[6 12 18]);

% T4
T4_con1=filtspikes(T4_noLEDspikes,0,'stimcond',[1 7 13]);
T4_con2=filtspikes(T4_noLEDspikes,0,'stimcond',[2 8 14]);
T4_con3=filtspikes(T4_noLEDspikes,0,'stimcond',[3 9 15]);
T4_con4=filtspikes(T4_noLEDspikes,0,'stimcond',[4 10 16]);
T4_con5=filtspikes(T4_noLEDspikes,0,'stimcond',[5 11 17]);
T4_con6=filtspikes(T4_noLEDspikes,0,'stimcond',[6 12 18]);

con1=[T1_con1.spiketimes T2_con1.spiketimes T3_con1.spiketimes T4_con1.spiketimes];
con2=[T1_con2.spiketimes T2_con2.spiketimes T3_con2.spiketimes T4_con2.spiketimes];
con3=[T1_con3.spiketimes T2_con3.spiketimes T3_con3.spiketimes T4_con3.spiketimes];
con4=[T1_con4.spiketimes T2_con4.spiketimes T3_con4.spiketimes T4_con4.spiketimes];
con5=[T1_con5.spiketimes T2_con5.spiketimes T3_con5.spiketimes T4_con5.spiketimes];
con6=[T1_con6.spiketimes T2_con6.spiketimes T3_con6.spiketimes T4_con6.spiketimes];

allCon1=concatSpikes_noAssigns(T1_con1,T2_con1);
allCon1=concatSpikes_noAssigns(allCon1,T3_con1);
allCon1=concatSpikes_noAssigns(allCon1,T4_con1);

allCon2=concatSpikes_noAssigns(T1_con2,T2_con2);
allCon2=concatSpikes_noAssigns(allCon2,T3_con2);
allCon2=concatSpikes_noAssigns(allCon2,T4_con2);

allCon3=concatSpikes_noAssigns(T1_con3,T2_con3);
allCon3=concatSpikes_noAssigns(allCon3,T3_con3);
allCon3=concatSpikes_noAssigns(allCon3,T4_con3);

allCon4=concatSpikes_noAssigns(T1_con4,T2_con4);
allCon4=concatSpikes_noAssigns(allCon4,T3_con4);
allCon4=concatSpikes_noAssigns(allCon4,T4_con4);

allCon5=concatSpikes_noAssigns(T1_con5,T2_con5);
allCon5=concatSpikes_noAssigns(allCon5,T3_con5);
allCon5=concatSpikes_noAssigns(allCon5,T4_con5);

allCon6=concatSpikes_noAssigns(T1_con6,T2_con6);
allCon6=concatSpikes_noAssigns(allCon6,T3_con6);
allCon6=concatSpikes_noAssigns(allCon6,T4_con6);

% psthWithOnlySpiketimes(con1,0.01,4,32000);
% h=axes;
% psth(allCon6,50,h,1,4);

h=axes;
[a,b,c,d,e,x6,y6]=psth(allCon6,100,h,1,4);
h=axes;
[a,b,c,d,e,x5,y5]=psth(allCon5,100,h,1,4);
h=axes;
[a,b,c,d,e,x4,y4]=psth(allCon4,100,h,1,4);
h=axes;
[a,b,c,d,e,x3,y3]=psth(allCon3,100,h,1,4);
h=axes;
[a,b,c,d,e,x2,y2]=psth(allCon2,100,h,1,4);
h=axes;
[a,b,c,d,e,x1,y1]=psth(allCon1,100,h,1,4);

baseInds=find(x1<1);
m1=mean(y1(baseInds));
m2=mean(y2(baseInds));
m3=mean(y3(baseInds));
m4=mean(y4(baseInds));
m5=mean(y5(baseInds));
m6=mean(y6(baseInds));
m=mean([m1 m2 m3 m4 m5 m6]);

figure();
plot(x1,y1-m1+m,'black');
hold on;
plot(x2,y2-m2+m,'green');
plot(x3,y3-m3+m,'blue');
plot(x4,y4-m4+m,'cyan');
plot(x5,y5-m5+m,'magenta');
plot(x6,y6-m6+m,'red');
% figure();
% plot(x1,mean([y1-m1+m y2-m2+m],2),'black');
% hold on;
% plot(x3,mean([y3-m3+m y4-m4+m],2),'blue');
% plot(x5,mean([y5-m5+m y6-m6+m],2),'red');