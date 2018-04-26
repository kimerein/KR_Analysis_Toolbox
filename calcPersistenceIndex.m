function [a,b,c,d,all_persbaseline35,all_persbaseline48]=calcPersistenceIndex(dataDir)

% dataDir='W:\Analysis Computer\Persistent Units\Mawake54 new units\T4\';

% a is evoked before led35, c is evoked persistent_window led35 (both minus
% led35 baseline
% b is evoked before led48, d is evoked persistent_window led48 (both minus
% led48 baseline 
files_before35=dir([dataDir 'T_*before35.mat']);
files_pers35=dir([dataDir 'T_*pers35.mat']);
files_persbaseline35=dir([dataDir 'T_*persbaseline35.mat']);
files_before=dir([dataDir 'T_*before.mat']);
files_pers=dir([dataDir 'T_*pers.mat']);
files_persbaseline48=dir([dataDir 'T_*persbaseline48.mat']);

FileList_before35={files_before35.name}';
FileList_pers35={files_pers35.name}';
FileList_persbaseline35={files_persbaseline35.name}';
FileList_before={files_before.name}';
FileList_pers={files_pers.name}';
FileList_persbaseline48={files_persbaseline48.name}';

disp(length(FileList_before35));
disp(length(FileList_pers35));
disp(length(FileList_persbaseline35));
disp(length(FileList_before));
disp(length(FileList_pers));
disp(length(FileList_persbaseline48));

all_before35=[];
all_pers35=[];
all_persbaseline35=[];
all_before=[];
all_pers=[];
all_persbaseline48=[];

for i=1:length(FileList_before35)
    a=load([dataDir FileList_before35{i}]);
    nb=a.nb;
    a=load([dataDir FileList_pers35{i}]);
    n=a.n;
    all_before35=[all_before35; mean(nb)];
    all_pers35=[all_pers35; mean(n)];
    
    a=load([dataDir FileList_persbaseline35{i}]);
    n=a.n;
    all_persbaseline35=[all_persbaseline35; mean(n)];
    a=load([dataDir FileList_persbaseline48{i}]);
    n=a.n;
    all_persbaseline48=[all_persbaseline48; mean(n)];
    
    a=load([dataDir FileList_before{i}]);
    nb=a.nb;
    a=load([dataDir FileList_pers{i}]);
    n=a.n;
    all_before=[all_before; mean(nb)];
    all_pers=[all_pers; mean(n)];
end

a=all_before35-all_persbaseline35;
c=all_pers35-all_persbaseline35;
b=all_before-all_persbaseline48;
d=all_pers-all_persbaseline48;

% Get before and pers corrcoeff
% files=dir([dataDir 'T_*before35.mat']);
% filesp=dir([dataDir 'T_*pers35.mat']);
% FileList={files.name}';
% FileListp={filesp.name}';
% % disp(FileList);
% cc_35=zeros(length(FileList),1);
% allnb35=[];
% alln35=[];
% for i=1:length(FileList)
%     a=load([dataDir FileList{i}]);
%     nb=a.nb;
%     a=load([dataDir FileListp{i}]);
%     n=a.n;
%     if addnoise==1
%         a=0;
%         b=0.0001;
%         nb=nb+(a + (b-a).*rand(size(nb)));
%         n=n+(a + (b-a).*rand(size(nb)));
%     end
%     disp(FileList{i});
%     disp(FileListp{i});
%     if shuffle==1
%         a=corrcoef(nb(randperm(length(nb))),n(randperm(length(n))));
%     else
%         a=corrcoef(nb,n);
%     end
%     cc_35(i)=a(1,2);
%     allnb35=[allnb35; nb];
%     alln35=[alln35; n];
% end    
% 
% cc_35(isnan(cc_35))=0;
% cc_pers(isnan(cc_pers))=0;
% 
% figure(); 
% scatter(cc_35,cc_pers);