function output=put_ArchT_units_together(filedir)

alignToDeepestUnit=1;
alignToBiggestChange=0;

listing=dir(filedir);

% Get folders
j=1;
for i=1:length(listing)
    if listing(i).isdir==1
        folders{j}=[filedir '\' listing(i).name];
        j=j+1;
    end
end

% Get files from folders
all_amp_change=[];
all_pvals=[];
all_visev_change=[];
all_colors=[];
all_depths=[];
all_halfWidths=[];
all_wvfms=[];
for i=1:length(folders)
    if exist([folders{i} '\output.mat'])
        a=load([folders{i} '\output.mat']);
        all_amp_change=[all_amp_change a.output.amp_change];
        curr_amp_change=a.output.amp_change;
        all_pvals=[all_pvals a.output.pvals];
        all_visev_change=[all_visev_change a.output.visev_change];
        all_colors=[all_colors; a.output.colors];
    end
    if exist([folders{i} '\all_width_and_depth.mat'])
        a=load([folders{i} '\all_width_and_depth.mat']);
        all_wvfms=[all_wvfms; a.all_wvfms];
        all_halfWidths=[all_halfWidths; a.all_halfWidths];
        
        if alignToDeepestUnit==1
%             a.all_depths=a.all_depths-max(a.all_depths(a.all_halfWidths>0.22*10^-3))+16;
            a.all_depths=a.all_depths-max(a.all_depths)+16;
%             figure();
%             scatter(a.all_halfWidths,a.all_depths);
        end
        if alignToBiggestChange==1
            [~,mi]=max(curr_amp_change);
            a.all_depths=a.all_depths-a.all_depths(mi)+16;
        end

        temp=folders{i};
        if strcmp(temp(end-4+1:end),'M355')
%             figure();
%             scatter(a.all_halfWidths,a.all_depths);
            a.all_depths=a.all_depths-16;
            a.all_depths=a.all_depths/2;
            a.all_depths=a.all_depths+16;
%             figure();
%             scatter(a.all_halfWidths,a.all_depths);
        end
        if strcmp(temp(end-4+1:end),'M374')
%             figure();
%             scatter(a.all_halfWidths,a.all_depths);
            a.all_depths=a.all_depths-16;
            a.all_depths=a.all_depths/2;
            a.all_depths=a.all_depths+16;
%             figure();
%             scatter(a.all_halfWidths,a.all_depths);
        end
        
        all_depths=[all_depths; a.all_depths];
    end
    
    
    
end

output.amp_change=all_amp_change;
output.pvals=all_pvals;
output.visev_change=all_visev_change;
output.colors=all_colors;
output.depths=all_depths;
output.halfWidths=all_halfWidths;
output.wvfms=all_wvfms;
plotResults(output);

end

function plotResults(output)

amp_change=output.amp_change;
visev_change=output.visev_change;
pvals=output.pvals;
cc=output.colors;
depths=output.depths;

cmap=colormap(jet(20));

figure();
for i=1:length(amp_change)
    if output.halfWidths(i)<=0.22*10^-3
%         scatter(amp_change(i),depths(i),[],cc(i,:),'filled');
    else
%     scatter(amp_change(i),visev_change(i),[],cc(i,:));
%     scatter(amp_change(i),pvals(i),[],cc(i,:));
%     scatter(amp_change(i),depths(i),[],cc(i,:));
%     scatter(visev_change(i),pvals(i),[],cc(i,:));
    scatter(amp_change(i),visev_change(i),[],cc(i,:));
    end
    hold on;
end
colorbar;
xlabel('amp change of suppression at led onset');
% ylabel('change in activity during visual stimulus');
ylabel('pval of change of suppression at led onset');

end