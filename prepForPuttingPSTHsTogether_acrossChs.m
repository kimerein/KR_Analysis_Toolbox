first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120221_Ch16_justEvoked.mat');
first1=first.y;
first2=first.y1;
first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120221_Ch16.mat');
first3=first.y2;
first4=first.y3;
first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120228_Ch16_justEvoked.mat');
second1=first.y;
second2=first.y1;
first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120228_Ch16.mat');
second3=first.y2;
second4=first.y3;
first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120330_Ch16.mat');
third1=first.y;
third2=first.y1;
third3=first.y2;
third4=first.y3;

% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120221_Ch2_justEvoked.mat');
% first5=first.y;
% first6=first.y1;
% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120221_Ch2.mat');
% first7=first.y2;
% first8=first.y3;
% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120228_Ch2_justEvoked.mat');
% second5=first.y;
% second6=first.y1;
% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120228_Ch2.mat');
% second7=first.y2;
% second8=first.y3;
% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120330_Ch2.mat');
% third5=first.y;
% third6=first.y1;
% third7=first.y2;
% third8=first.y3;

% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120221_Ch4_justEvoked.mat');
% first9=first.y;
% first10=first.y1;
% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120221_Ch4.mat');
% first11=first.y2;
% first12=first.y3;
% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120228_Ch4_justEvoked.mat');
% second9=first.y;
% second10=first.y1;
% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120228_Ch4.mat');
% second11=first.y2;
% second12=first.y3;
% first=load('C:\Users\Kim\Documents\MATLAB\Spikes Structs\Layer-Specific PSTHs\120330_Ch4.mat');
% third9=first.y;
% third10=first.y1;
% third11=first.y2;
% third12=first.y3;

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
% dataY1(4,:)=first5;
% dataY1(5,:)=second5;
% dataY1(6,:)=third5;
% dataY1(7,:)=first9;
% dataY1(8,:)=second9;
% dataY1(9,:)=third9;

dataY2(1,:)=first2;
dataY2(2,:)=second2;
dataY2(3,:)=third2;
% dataY2(4,:)=first6;
% dataY2(5,:)=second6;
% dataY2(6,:)=third6;
% dataY2(7,:)=first10;
% dataY2(8,:)=second10;
% dataY2(9,:)=third10;
% putPSTHsTogetherFromAcrossFiles(second.x,dataY1,dataY2,30,0.05,0.001);
putPSTHsTogetherFromAcrossFiles(second.x,dataY1,dataY2,30,0.03,0.001);

clear dataY1 dataY2;

dataY1(1,:)=first3;
dataY1(2,:)=second3;
dataY1(3,:)=third3;
dataY2(1,:)=first4;
dataY2(2,:)=second4;
dataY2(3,:)=third4;
% putPSTHsTogetherFromAcrossFiles(second.x,dataY1,dataY2,30,0.05,0.001);
% putTausTogetherFromAcrossFiles(second.x,dataY1,dataY2,10,0.05,0.001);