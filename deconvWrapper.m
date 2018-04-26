function deconvWrapper(headerDir,newDir)

listing=dir(headerDir);

for i=3:length(listing)
    currFolder=listing(i).name;
    if ~exist([headerDir '\' currFolder '\psths.mat'],'file')
    else
        a=load([headerDir '\' currFolder '\psths.mat']);
        curr=a.psths;
        deconvolutionDynamics(currFolder,[],[],curr,newDir);
    end
end

