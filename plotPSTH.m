function [fig_handles,axes_handles,colorMapCode]=plotPSTH(spikes,params,ONresponseWindow,OFFresponseWindow,quantResponse)
% ONperiod is a 2-element vector indicating the beginning and end of the
% region to consider for quantifying the ON response
% OFFperiod is a 2-element vector indicating the beginning and end of the
% region to consider for quantifying the OFF response
% quantResponse is the method for quantifying the ON or OFF response
%       can be 'peak','integral' or 'average'
%       Method - Response is the ... 
%       'peak' - peak firing rate during the specified period
%       'integral' - area under the firing rate curve of the specified
%       period
%       'average' - average of the firing rate across the entire specified
%       period
% fig_handles is a structure array containing handles to the different PSTH
% figures, with associated figure descriptions

fig_handles=[];
axes_handles=[];

f1=figure;
fig_handles=[fig_handles; f1];
% Plan figure layout
nfigRows=1;
nfigCols=3;
figParams=[];
figParams.matpos=[0.05 0.05 0.05 0.05];
figParams.figmargin=[0.05 0.05 0.05 0.05];
figParams.matmargin=[0.05 0.05 0.05 0.05];
figParams.cellmargin=[0.1 0.1 0.1 0.1];

% Assign different colors to different stimulus conditions
[lineColors,colorMapCode]=makeColorCode(params);

% Plot the total PSTH, collapsed over all stimulus conditions
curr_h=axesmatrix(nfigRows,nfigCols,1);
%[curr_h,maxFR,binWidth,ONresponse,OFFresponse]=superPSTH(hax, spikes, 'set', 0.05, lineColors, params, 'Total PSTH - All Spikes, All Stimulus Conditions',[],quantResponse,ONresponseWindow,OFFresponseWindow);
%K_psth(spikes,0.5,params.totalTrialLength,'k');
binWidth=0.05;   % binWidth in s
%K_psth(spikes,binWidth,params.totalTrialLength,'k'); 
[h,maxFR,binWidth,ONresponse,OFFresponse]=K_superPSTH(spikes,'set',0.115,[0 0 0],0,params,'Total PSTH',[],quantResponse,ONresponseWindow,OFFresponseWindow);
s=struct('figure_handle',curr_h,'figure_description','total PSTH, collapsed over all stimulus conditions');
axes_handles=[axes_handles; s];
%setaxesOnaxesmatrix(hax,nfigRows,nfigCols,1,figParams,f1);
%setaxesOnaxesmatrix(curr_h,nfigRows,nfigCols,1,[],curr_h);

% Plot PSTH collapsed over Stimulus Variable 2
%set(f1,'CurrentAxes',axesmatrix(nfigRows,nfigCols,2,figParams));
curr_h=axesmatrix(nfigRows,nfigCols,2);
someSpikes=[];
for i=1:length(params.Var1_values) 
    someSpikes=[someSpikes; filtspikes(spikes,0,'stimcond',i:1:i+length(params.Var2_values)-1)];
end
[h,maxFR,binWidth,ONresponse,OFFresponse]=K_superPSTH(someSpikes,'set',0.115,lineColors,length(params.Var1_values)*length(params.Var2_values),params,['PSTH - Effects of ' params.Var1_name],[],quantResponse,ONresponseWindow,OFFresponseWindow);
%[curr_h,maxFR,binWidth,ONresponse,OFFresponse]=superPSTH(hax,someSpikes, 'set', 0.05, lineColors, params, ['PSTH - Effects of ' params.Var1_name],[],quantResponse,ONresponseWindow,OFFresponseWindow);
s=struct('figure_handle',curr_h,'figure_description',['Superimposed: effects of ' params.Var1_name ' - PSTH collapsed over ' params.Var2_name]);
axes_handles=[axes_handles; s];

% Plot PSTH collapsed over Stimulus Variable 1
%set(f1,'CurrentAxes',axesmatrix(nfigRows,nfigCols,3,figParams,f1));
curr_h=axesmatrix(nfigRows,nfigCols,3);
someSpikes=[];
for i=1:length(params.Var2_values) 
    someSpikes=[someSpikes; filtspikes(spikes,0,'stimcond',i:length(params.Var2_values):length(params.Var1_values)*length(params.Var2_values))];
end
[h,maxFR2,binWidth,ONresponse,OFFresponse]=K_superPSTH(someSpikes,'set',0.115,lineColors,length(params.Var1_values)*length(params.Var2_values)+length(params.Var1_values),params,['PSTH - Effects of ' params.Var2_name],[],quantResponse,ONresponseWindow,OFFresponseWindow);
if maxFR2>maxFR
    maxFR=maxFR2;
end
s=struct('figure_handle',curr_h,'figure_description',['Superimposed: effects of ' params.Var2_name ' - PSTH collapsed over ' params.Var1_name]);
axes_handles=[axes_handles; s];

f2=figure;
fig_handles=[fig_handles; f2];
% Plan figure layout
nfigRows=1+length(params.Var1_values);
nfigCols=2+length(params.Var2_values);
figParams.matpos=[0.05 0.05 0.05 0.05];
figParams.figmargin=[0.05 0.05 0.05 0.05];
figParams.matmargin=[0.05 0.05 0.05 0.05];
figParams.cellmargin=[0.1 0.1 0.1 0.1];

% For each row of PSTH figures, at end of row, show ON and OFF responses
% in bar graph, with order of bars corresponding to the order of graphs in
% that row
ONresponseArray=zeros(nfigRows,nfigCols);
OFFresponseArray=zeros(nfigRows,nfigCols);

k=2;
% Collapsed over Var1
for i=1:length(params.Var2_values) 
    %set(f2,'CurrentAxes',axesmatrix(nfigRows,nfigCols,k,figParams,f2));
    curr_h=axesmatrix(nfigRows,nfigCols,k);
    someSpikes=filtspikes(spikes,0,'stimcond',i:length(params.Var2_values):length(params.Var1_values)*length(params.Var2_values));
    [h,maxFR2,binWidth,ONresponse,OFFresponse]=K_superPSTH(someSpikes,'set',binWidth,lineColors,length(params.Var1_values)*length(params.Var2_values)+length(params.Var1_values),params,['Efx. of ' params.Var2_name],[0 params.totalTrialLength 0 maxFR],quantResponse,ONresponseWindow,OFFresponseWindow);
    s=struct('figure_handle',curr_h,'figure_description',['effects of ' params.Var2_name ' - PSTH collapsed over ' params.Var1_name]);
    axes_handles=[axes_handles; s];
    ro=floor(k/nfigCols)+1;
    co=mod(k,nfigCols);
    if co==0
        co=nfigCols;
    end 
    ONresponseArray(ro,co)=ONresponse(1);
    OFFresponseArray(ro,co)=OFFresponse(1);
    k=k+1;
end
orientONs=zeros(length(params.Var1_values),1);
orientOFFs=zeros(length(params.Var1_values),1);
% Collapsed over Var2
for i=1:length(params.Var1_values)
    %set(f2,'CurrentAxes',axesmatrix(nfigRows,nfigCols,i*nfigCols,figParams,f2));
    curr_h=axesmatrix(nfigRows,nfigCols,i*nfigCols+1);
    someSpikes=filtspikes(spikes,0,'stimcond',i:1:i+length(params.Var2_values)-1);
    [h,maxFR2,binWidth,ONresponse,OFFresponse]=K_superPSTH(someSpikes,'set',binWidth,lineColors,length(params.Var1_values)*length(params.Var2_values),params,['Efx. of ' params.Var1_name],[0 params.totalTrialLength 0 maxFR],quantResponse,ONresponseWindow,OFFresponseWindow);
    s=struct('figure_handle',curr_h,'figure_description',['effects of ' params.Var1_name ' - PSTH collapsed over ' params.Var2_name]);
    axes_handles=[axes_handles; s];
    orientONs(i)=ONresponse;
    orientOFFs(i)=OFFresponse;
    ro=floor((i*nfigCols+1)/nfigCols)+1;
    co=mod(i*nfigCols+1,nfigCols);
    if co==0
        co=nfigCols;
    end
    ONresponseArray(ro,co)=ONresponse(1);
    OFFresponseArray(ro,co)=OFFresponse(1);
end
% Separate plots for each stimulus condition (combo. of Var1 and Var2)
k=nfigCols+2;
keepTheseAxes=[];
for i=1:length(params.Var1_values)
    for j=1:length(params.Var2_values)
        %set(f2,'CurrentAxes',axesmatrix(nfigRows,nfigCols,(i-1)*length(params.Var2_values)+j,figParams,f2));
        curr_h=axesmatrix(nfigRows,nfigCols,k);
        keepTheseAxes=[keepTheseAxes curr_h];
        someSpikes=filtspikes(spikes,0,'stimcond',(i-1)*length(params.Var2_values)+j);
        [h,maxFR2,binWidth,ONresponse,OFFresponse]=K_superPSTH(someSpikes,'set',binWidth,lineColors,0,params,[num2str(params.Var1_values(i)) ' ' num2str(params.Var2_values(j))],[0 params.totalTrialLength 0 maxFR],quantResponse,ONresponseWindow,OFFresponseWindow);
        if maxFR2>maxFR
            maxFR=maxFR2;
        end
        s=struct('figure_handle',curr_h,'figure_description',[params.Var1_name ' is ' num2str(params.Var1_values(i)) ', ' params.Var2_name ' is' num2str(params.Var2_values(j))]);
        axes_handles=[axes_handles; s];
        ro=floor(k/nfigCols)+1;
        co=mod(k,nfigCols);
        if co==0
            co=nfigCols;
        end
        ONresponseArray(ro,co)=ONresponse(1);
        OFFresponseArray(ro,co)=OFFresponse(1);
        if mod(k,nfigCols)==nfigCols-1
            k=k+3;
        else
            k=k+1;
        end 
    end
end
for i=1:length(keepTheseAxes)
    axis(keepTheseAxes(i),[0 params.totalTrialLength 0 maxFR]);
end

% Plot ON and OFF responses
for i=1:nfigRows
    %set(f2,'CurrentAxes',axesmatrix(nfigRows,nfigCols,nfigCols*i,figParams,f2));
    curr_h=axesmatrix(nfigRows,nfigCols,nfigCols*i);
    h=bar(curr_h,[ONresponseArray(i,:); OFFresponseArray(i,:)],1,'group');
    if any(any(ONresponseArray~=0,1),2) || any(any(OFFresponseArray~=0,1),2)
        ylim([0 max([ONresponseArray(i,:) OFFresponseArray(i,:)])]);
    end
    title('ON,OFF');
end

% If stimulus variable 1 is Orientation and stimulus variable 2 is Phase,
% make polar plots for ON and OFF responses
f3=figure;
fig_handles=[fig_handles; f3];
if strcmp(params.Var1_name,'Orientation') && ~isempty(orientONs)
    orientStep=360/length(params.Var1_values);
    orients=0:orientStep:360;
    orientONs=[orientONs; orientONs(1)];
    orientOFFs=[orientOFFs; orientOFFs(1)];
    ONnormMags=orientONs/max(orientONs);
    OFFnormMags=orientOFFs/max(orientOFFs);
    curr_h=subplot(1,2,1);
    s=struct('figure_handle',curr_h,'figure_description','Orientation tuning - polar plot');
    axes_handles=[axes_handles; s];
    if max(orientONs)>max(orientOFFs)
        polar(deg2rad(orients)',orientONs,'-g');
        hold on;
        title('Orientation Tuning');
        polar(deg2rad(orients)',orientOFFs,'-r');
        curr_h=subplot(1,2,2);
        s=struct('figure_handle',curr_h,'figure_description','Orientation tuning - polar plot');
        axes_handles=[axes_handles; s];
        polar(deg2rad(orients)',ONnormMags,'-g');
        hold on;
        polar(deg2rad(orients)',OFFnormMags,'-r');
        title('ON: Green, OFF: Red');
    else
        polar(deg2rad(orients)',orientOFFs,'-r');
        hold on;
        title('Orientation Tuning');
        polar(deg2rad(orients)',orientONs,'-g');
        curr_h=subplot(1,2,2);
        s=struct('figure_handle',curr_h,'figure_description','Orientation tuning - Normalized polar plot');
        axes_handles=[axes_handles; s];
        polar(deg2rad(orients)',OFFnormMags,'-r');
        hold on;
        polar(deg2rad(orients)',ONnormMags,'-g');
        title('ON: Green, OFF: Red');
    end
end


