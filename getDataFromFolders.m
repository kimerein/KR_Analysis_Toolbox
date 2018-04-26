function allData=getDataFromFolders(headerDir,dataName,inFileName,vertcat)

listing=dir(headerDir);

allData=[];
for i=3:length(listing)
    currFolder=listing(i).name;
    if ~exist([headerDir '\' currFolder '\' dataName '.mat'],'file')
    else
        a=load([headerDir '\' currFolder '\' dataName '.mat']);
        curr=a.(inFileName);
        if vertcat==1
            allData=[allData; curr];
        else
            allData=[allData curr];
        end
    end
end