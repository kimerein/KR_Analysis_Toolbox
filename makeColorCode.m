function [lineColors,colorMapCode]=makeColorCode(params)
% Create a color map for different stimulus conditions
% Different color for each combination of Var1 and Var2 values
% PLUS different colors for data collapsed over Var1 or Var2
% Data collapsed across Var1 AND Var2 is plotted in black
%
% params are the stimulus condition parameters
% params is a structure that must contain
%   params.Var1_values = a vector of stimulus parameter 1 values
%   params.Var2_values = a vector of stimulus parameter 2 values
% 
% colorMapCode
% This structure array contains three fields:
% color = three-element color row vector with values between 0 and 1
% var1_values = vector of stimulus parameter 1 values (in conjunction with 
%   the var2_values) represented by this color
% var2_values = vector of stimulus parameter 2 values (in conjunction with 
%   the var1_values) represented by this color
%
% lineColors
% an n-by-3 matrix where each row is an RGB color value (each
% element between 0 and 1)
% the number of rows is the number of stimulus conditions + the number of
% conditions collapsed over stim. condition varying parameter 1 + the
% number of conditions collapsed over stim. condition varying parameter 2
% Essentially, this matrix specifies a color for each condition to be
% plotted in plotLFPbyStimCond.m

nConditions=length(params.Var1_values)*length(params.Var2_values)+length(params.Var1_values)+length(params.Var2_values);
cIt=ceil(nConditions^(1/3));
m=0:0.8/(cIt-1):0.8;
j=0:0.8/(cIt-1):0.8;
k=0:0.8/(cIt-1):0.8;
lineColors=zeros(length(m)^3,3);
l=1;
for m1=m
    for j1=j
        for k1=k
            lineColors(l,:)=[m1 j1 k1];
            l=l+1;
        end
    end
end
lineColors=lineColors(randperm(size(lineColors,1)),:);
% Create a color map code = colorMapCode
% This structure array contains three fields:
% color = three-element color row vector with values between 0 and 1
% var1_values = vector of stimulus parameter 1 values (in conjunction with 
%   the var2_values) represented by this color
% var2_values = vector of stimulus parameter 2 values (in conjunction with 
%   the var1_values) represented by this color
colorMapCode=[];
l=1;
for i=1:length(params.Var1_values)
    for j=1:length(params.Var2_values)
        c=struct('color',lineColors(l,:),'var1_values',params.Var1_values(i),'var2_values',params.Var2_values(j));
        colorMapCode=[colorMapCode; c];
        l=l+1;
    end
end
for i=1:length(params.Var1_values)
    c=struct('color',lineColors(l,:),'var1_values',params.Var1_values(i),'var2_values',params.Var2_values(:));
    colorMapCode=[colorMapCode; c];  
    l=l+1;
end
for i=1:length(params.Var2_values)
    c=struct('color',lineColors(l,:),'var1_values',params.Var1_values(:),'var2_values',params.Var2_values(i));
    colorMapCode=[colorMapCode; c];  
    l=l+1;
end