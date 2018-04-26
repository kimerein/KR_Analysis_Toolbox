function concatWidthAndDepth(filedir,trode_numbers)

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
all_wvfms=[];
all_halfWidths=[];
all_depths=[];
all_assigns=[];
if length(folders(3:end))~=length(trode_numbers)
    disp('trode numbers length does not match number of folders in this directory');
    return
end
j=1;
for i=1:length(folders)
    if exist([folders{i} '\width_and_depth.mat'])
        a=load([folders{i} '\width_and_depth.mat']);
        all_wvfms=[all_wvfms; a.unit_wvfms];
        all_halfWidths=[all_halfWidths; a.unit_halfWidths];
        all_depths=[all_depths; a.unit_depths+(trode_numbers(j)-1)*4];
        all_assigns=[all_assigns a.useAssigns];
        j=j+1;
    end
end

save([filedir '\all_width_and_depth.mat'],'all_wvfms','all_halfWidths','all_depths','all_assigns');

