function plotNtsr1Scatter(datadir)

takeFSInsteadOfRS=false;

if iscell(datadir)
    ampsAcrossCells=[];
    pvalsAcrossCells=[];
    exptGrpAcrossCells=[];
    isNtsr1AcrossCells=[];
    ZdepthAcrossCells=[];
    allUsedWvfms=[];
    for i=1:length(datadir)
        d=datadir{i};
        a=load([d '\' 'output.mat']);
        ntsr=a.output;
        a=load([d '\' 'wvfms.mat']);
        wvfms=a.wvfms;
        a=load([d '\' 'classifyAsNtsr1.mat']);
        classifyAsNtsr1=a.classifyAsNtsr1;
        a=load([d '\' 'deepenough.mat']);
        deepenough=a.deepenough;
        if length(deepenough)==1
            deepenough=ones(1,length(wvfms.hw)).*deepenough;
        end
        a=load([d '\' 'exptType.mat']);
        exptType=a.exptType;
        % topCh_wrt_bottomMostChWSpikes is in reverse order, from deepest
        % to most superficial; ch 1 is the bottom-most ch w spikes across
        % the entire array
        a=load([d '\' 'topCh_wrt_bottomMostChWSpikes.mat']);
        topCh_wrt_bottomMostChWSpikes=a.topCh_wrt_bottomMostChWSpikes;
        a=load([d '\' 'ch_spacing.mat']);
        ch_spacing=a.ch_spacing;
        if exist([d '\' 'unit_depths.mat'],'file')
            alreadyExists=true;
            a=load([d '\' 'unit_depths.mat']);
            unit_depths=a.unit_depths;
        else
            alreadyExists=false;
            unit_depths=[];
        end
        if isempty(unit_depths)
            unit_depths=getUnitDepths(wvfms);
        end
        if alreadyExists==false
            save([d '\' 'unit_depths.mat'],'unit_depths');
        end
        % calculate real unit depth
        distance_real_Z=calcRealUnitDepth(unit_depths,topCh_wrt_bottomMostChWSpikes,ch_spacing);
        % take unit type
        a=load([d '\' 'isFs.mat']);
        isFs=a.isFs;
        if takeFSInsteadOfRS==true
            takeWvfms=wvfms.spikewvfms(isFs==1,:,:);
            classifyAsNtsr1=classifyAsNtsr1(isFs==1);
            distance_real_Z=distance_real_Z(isFs==1);
            deepenough=deepenough(isFs==1);
            f=fieldnames(ntsr);
            for j=1:length(f)
                temp=ntsr.(f{j});
                if length(temp)~=length(isFs)
                    continue
                end
                temp=temp(isFs==1);
                ntsr.(f{j})=temp;
            end
        else
%             takeWvfms=wvfms.spikewvfms(isFs==0,:,:);
            takeWvfms=wvfms.spikewvfms(isFs==0 & classifyAsNtsr1==1,:,:);
            classifyAsNtsr1=classifyAsNtsr1(isFs==0);
            distance_real_Z=distance_real_Z(isFs==0);
            deepenough=deepenough(isFs==0);
            f=fieldnames(ntsr);
            for j=1:length(f)
                temp=ntsr.(f{j});
                if length(temp)~=length(isFs)
                    continue
                end
                temp=temp(isFs==0);
                ntsr.(f{j})=temp;
            end
        end        
        % get what to plot
        switch exptType
            case 'ChR'
                plotThisAmp='quick_amp_change';
                plotThisPval='quick_change_pvals';
                exptGrpAcrossCells=[exptGrpAcrossCells ones(size(classifyAsNtsr1)).*1];
            case 'ArchT'
                if isfield(ntsr,'quick_amp_change')
                    plotThisAmp='quick_amp_change';
                    plotThisPval='quick_change_pvals';
                else
                    plotThisAmp='amp_change';
                    plotThisPval='pvals';
                end
                exptGrpAcrossCells=[exptGrpAcrossCells ones(size(classifyAsNtsr1)).*2];
            case 'Gtacr2'
                plotThisAmp='quick_amp_change';
                plotThisPval='quick_change_pvals';
                exptGrpAcrossCells=[exptGrpAcrossCells ones(size(classifyAsNtsr1)).*3];
            otherwise 
                error('unrecognized exptType');
        end
        ampsAcrossCells=[ampsAcrossCells ntsr.(plotThisAmp)];
        pvalsAcrossCells=[pvalsAcrossCells ntsr.(plotThisPval)];
        isNtsr1AcrossCells=[isNtsr1AcrossCells classifyAsNtsr1];
        ZdepthAcrossCells=[ZdepthAcrossCells distance_real_Z];
        if length(size(takeWvfms))<3
        else
            [~,mi]=nanmin(nanmin(takeWvfms,[],2),[],3);
            if isempty(allUsedWvfms) || size(takeWvfms,2)==size(allUsedWvfms,2)
                for k=1:size(takeWvfms,1)
                    allUsedWvfms=[allUsedWvfms; reshape(takeWvfms(k,:,mi(k)),1,size(takeWvfms,2))];
                end
            end
%             if (size(allUsedWvfms,2)==0) || size(takeWvfms,2)*size(takeWvfms,3)==size(allUsedWvfms,2)
%                 allUsedWvfms=[allUsedWvfms; reshape(takeWvfms,size(takeWvfms,1),size(takeWvfms,2)*size(takeWvfms,3))];
%             end
        end
    end
else
    error('expected datadir to be a cell array');
end

% Plot wvfms figure
plotWStderr(0:1/25000:(size(allUsedWvfms,2)-1)*(1/25000),allUsedWvfms,[],'k',[],0);

% Plot figure
sz=50;
figure();
for i=1:length(ampsAcrossCells)
    if exptGrpAcrossCells(i)==1
        c=[0 0.5 0.8];
    elseif exptGrpAcrossCells(i)==2
        c='r';
    elseif exptGrpAcrossCells(i)==3
        c='r';
    end
    if isNtsr1AcrossCells(i)==1
        m='k';
    else
        m=c;
    end
    if pvalsAcrossCells(i)<=0.05
        scatter(ampsAcrossCells(i),ZdepthAcrossCells(i),sz+10,c,'filled','MarkerEdgeColor',m,'LineWidth',1.5);
    else
        scatter(ampsAcrossCells(i),ZdepthAcrossCells(i),sz,c,'LineWidth',1.5);
    end
    hold on;
end

end

function plotWStderr(x,y1,y2,c1,c2,doFill)

figure();
plot(x,nanmean(y1,1),'Color',c1);
hold on;
if doFill==1
    fill([x fliplr(x)],[nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)) fliplr(nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)))],[0.5 0.5 0.5]);
end
plot(x,nanmean(y1,1),'Color',c1);
plot(x,nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);
plot(x,nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);

if ~isempty(y2)
    plot(x,nanmean(y2,1),'Color',c2);
    if doFill==1
        fill([x fliplr(x)],[nanmean(y2,1)+nanstd(y2,[],1)./sqrt(size(y2,1)) fliplr(nanmean(y2,1)-nanstd(y2,[],1)./sqrt(size(y2,1)))],[0.1 0.7 0.5]);
    end
    plot(x,nanmean(y2,1),'Color',c2);
    plot(x,nanmean(y2,1)+nanstd(y2,[],1)./sqrt(size(y2,1)),'Color',c2);
    plot(x,nanmean(y2,1)-nanstd(y2,[],1)./sqrt(size(y2,1)),'Color',c2);
    
    plot(x,nanmean(y1,1),'Color',c1);
end
    
end

function distance_real_Z=calcRealUnitDepth(unit_depths,topCh_wrt_bottomMostChWSpikes,ch_spacing)

% probe angle 40 deg from tangent to pia
assumeProbeAngle=40; % in degrees
% thus real z depth into cortex is depth on probe times sind(40)

% first get channel depth on probe wrt bottom-most channel with spikes
ch_depth_on_probe=topCh_wrt_bottomMostChWSpikes-unit_depths+1; % smallest number is deepest ch

% then get real distance along probe
distance_on_probe=(ch_depth_on_probe-1)*ch_spacing; % distance from bottom-most spike in microns

% take into account probe angle
distance_real_Z=distance_on_probe*sind(assumeProbeAngle); % microns from bottom-most spike

if size(distance_real_Z,1)>1
    distance_real_Z=distance_real_Z';
end

end

function unit_depths=getUnitDepths(wvfms)

unit_depths=nan(1,size(wvfms.spikewvfms,1));
for i=1:size(wvfms.spikewvfms,1)
    ampPerCh=nan(1,size(wvfms.spikewvfms,3));
    for j=1:size(wvfms.spikewvfms,3)
        ampPerCh(j)=nanmin(wvfms.spikewvfms(i,:,j));
    end
    ampPerCh=abs(ampPerCh);
    weightedCh=0;
    for j=1:size(wvfms.spikewvfms,3)
        weightedCh=weightedCh+j*ampPerCh(j)/(nansum(ampPerCh));
    end
    unit_depths(i)=weightedCh;
end

end