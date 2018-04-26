function [cc_pers,cc_35,allnb,alln,allnb35,alln35]=compareUnitCorrcoeffs()

shuffle=1;
addnoise=1;
makeProbFR=1;

dataDir='W:\Analysis Computer\Persistent Units\Mawake54 new units\T4\';

% Get before and pers corrcoeff
files=dir([dataDir 'T_*before.mat']);
filesp=dir([dataDir 'T_*pers.mat']);
FileList={files.name}';
FileListp={filesp.name}';
% disp(FileList);
cc_pers=zeros(length(FileList),1);
allnb=[];
alln=[];
for i=1:length(FileList)
    a=load([dataDir FileList{i}]);
    nb=a.nb;
    a=load([dataDir FileListp{i}]);
    n=a.n;
    if addnoise==1
        a=0;
        b=0.0001;
        nb=nb+(a + (b-a).*rand(size(nb)));
        n=n+(a + (b-a).*rand(size(nb)));
    end
    disp(FileList{i});
    disp(FileListp{i});
    if shuffle==1
        a=corrcoef(nb(randperm(length(nb))),n(randperm(length(n))));
    else
        a=corrcoef(nb,n);
    end
    cc_pers(i)=a(1,2);
    allnb=[allnb; nb];
    alln=[alln; n];
end

% Get before and pers corrcoeff
files=dir([dataDir 'T_*before35.mat']);
filesp=dir([dataDir 'T_*pers35.mat']);
FileList={files.name}';
FileListp={filesp.name}';
% disp(FileList);
cc_35=zeros(length(FileList),1);
allnb35=[];
alln35=[];
for i=1:length(FileList)
    a=load([dataDir FileList{i}]);
    nb=a.nb;
    a=load([dataDir FileListp{i}]);
    n=a.n;
    if addnoise==1
        a=0;
        b=0.0001;
        nb=nb+(a + (b-a).*rand(size(nb)));
        n=n+(a + (b-a).*rand(size(nb)));
    end
    disp(FileList{i});
    disp(FileListp{i});
    if shuffle==1
        a=corrcoef(nb(randperm(length(nb))),n(randperm(length(n))));
    else
        a=corrcoef(nb,n);
    end
    cc_35(i)=a(1,2);
    allnb35=[allnb35; nb];
    alln35=[alln35; n];
end    

cc_35(isnan(cc_35))=0;
cc_pers(isnan(cc_pers))=0;

figure(); 
scatter(cc_35,cc_pers);