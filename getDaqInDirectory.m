function [theseAreDaqs,LED,Trig]=getDaqInDirectory(expt_dir)

sortNames=1;
getLEDCond=1;
getTrigCond=1;
LED=[];
Trig=[];

ls=dir(expt_dir);
k=1;
theseAreDaqs=[];
for i=1:length(ls)
    thisname=ls(i).name;
    if ~isempty(regexp(thisname,'\.daq'))
        % is daq file
        % return name
        theseAreDaqs{k}=thisname;
        k=k+1;
    end
end

% sort names by number before daq
if sortNames==1
    daqNum=nan(1,length(theseAreDaqs));
    for i=1:length(theseAreDaqs)
        thisname=theseAreDaqs{i};
        r1=regexp(thisname,'_');
        if length(r1)>1
            r1=r1(2); % only for arbora's data
        end
        r2=regexp(thisname,'\.');
        daqNum(i)=str2num(thisname(r1+1:r2-1));
    end
    [~,sorti]=sort(daqNum);
    theseAreDaqs=theseAreDaqs(sorti);    
end

if getLEDCond==1
    LED=getCond(theseAreDaqs,'_LEDCond',expt_dir,'LEDCond');
end
if getTrigCond==1
    Trig=getCond(theseAreDaqs,'_TrigCond',expt_dir,'CondNum');    
end

end

function out=getCond(daqNames,fileNameEnd,dirName,nameInFile)

out=[];
for i=1:length(daqNames)
    thisname=daqNames{i};
    r=regexp(thisname,'\.');
    a=load([dirName '\' thisname(1:r-1) fileNameEnd '.mat']);
    out=[out a.(nameInFile)];
end

end
