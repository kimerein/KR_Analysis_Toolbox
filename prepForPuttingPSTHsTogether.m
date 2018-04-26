first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120221_Ch14_justEvoked.mat');
first1=first.y;
first2=first.y1;
first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120221_Ch14.mat');
first3=first.y2;
first4=first.y3;
first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120228_Ch14_justEvoked.mat');
second1=first.y;
second2=first.y1;
first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120228_Ch14.mat');
second3=first.y2;
second4=first.y3;
first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120330_Ch14.mat');
third1=first.y;
third2=first.y1;
third3=first.y2;
third4=first.y3;

% first1=first.y;
% first2=first.y1;
% first3=first.y2;
% first4=first.y3;
% 
% second1=second.y;
% second2=second.y1;
% second3=second.y2;
% second4=second.y3;
% 
% third1=third.y;
% third2=third.y1;
% third3=third.y2;
% third4=third.y3;

clear dataY1 dataY2;

dataY1(1,:)=first1;
dataY1(2,:)=second1;
dataY1(3,:)=third1;
dataY2(1,:)=first2;
dataY2(2,:)=second2;
dataY2(3,:)=third2;
% putPSTHsTogetherFromAcrossFiles(second.x,dataY1,dataY2,30,0.05,0.001);
putTausTogetherFromAcrossFiles(second.x,dataY1,dataY2,5,0.03,0.001);

clear dataY1 dataY2;

dataY1(1,:)=first3;
dataY1(2,:)=second3;
dataY1(3,:)=third3;
dataY2(1,:)=first4;
dataY2(2,:)=second4;
dataY2(3,:)=third4;
% putPSTHsTogetherFromAcrossFiles(second.x,dataY1,dataY2,30,0.05,0.001);
% putTausTogetherFromAcrossFiles(second.x,dataY1,dataY2,10,0.05,0.001);