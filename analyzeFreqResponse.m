function [F1,F2,F3,F2MinusF3,thalResults]=analyzeFreqResponse(useDir,listing)

putTogetherThalData=0;
putTogetherV1Data=1;
subtractNonspec=0;

doAmp=1;
doMU=0;
doUnitsOneByOne=1;
integralAlign=0;
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

allUnits=[];
fromMouse=[];
F1=[];
F2=[];
F3=[];
F2MinusF3=[];
thalResults=[];
for i=1:length(listing)
    currFolder=listing(i).name;
    if doMU==1
        if ~exist([useDir '\' currFolder '\MUcrosscorr_1.mat'],'file')
            allSU=[];
        else
            a=load([useDir '\' currFolder '\MUcrosscorr_1.mat']);
            allSU=a.MUcrosscorr;
        end
        if ~exist([useDir '\' currFolder '\responseFrequencies.mat'],'file')
            bigFres=[];
        else
            a=load([useDir '\' currFolder '\responseFrequencies.mat']);
            bigFres=a.bigFres;
        end
    else
        if ~exist([useDir '\' currFolder '\allSU1.mat'],'file')
            allSU=[];
        else
            a=load([useDir '\' currFolder '\allSU1.mat']);
            allSU=a.allSU;
            if putTogetherThalData==1
                fromMouse=[fromMouse i*ones(1,length(allSU))];
            end
        end
        if ~exist([useDir '\' currFolder '\bigFres1.mat'],'file')
            bigFres=[];
        else
            a=load([useDir '\' currFolder '\bigFres1.mat']);
            bigFres=a.bigFres;
        end
    end
    allUnits=[allUnits allSU];
end

if doAmp==1
    for i=1:length(allUnits)
        allUnits(i).responses{1}=sqrt(allUnits(i).responses{1});
    end
end

if putTogetherV1Data==1
    bigFres=bigFres(1,:);
    % Peak-norm. units
    avResp=zeros(size(allUnits(1).responses{1}));
    for i=1:length(allUnits)
        allUnits(i).responses{1}=allUnits(i).responses{1}./max(max(allUnits(i).responses{1}));
        avResp=avResp+allUnits(i).responses{1};
    end
    figure(); 
    imagesc(avResp(:,3:end).*2);
    % Get F1,F2,F3
    m=zeros(15,122);
    for i=1:length(allUnits)
        [F1(i,:),F2(i,:),F3(i,:),~,mcurr]=plotFreqResponseWithHarms(allUnits(i).responses{1},bigFres,0,0,subtractNonspec);
        F2MinusF3(i,:)=F2(i,:)-F3(i,:);
        m=m+mcurr;
    end
    figure(); imagesc(m);
    % Plot results
    figure();
    hax=axes();
    hl=semilogx(freqs,nanmean(F1,1),'Color','k');
    hold on;
    addErrBar(freqs,nanmean(F1,1),nanstd(F1,[],1)./sqrt(size(F1,1)),'y',hax,hl);
    hl=semilogx(freqs,nanmean(F2MinusF3,1),'Color',[0.8 0.8 0.8]);
    hold on;
    addErrBar(freqs,nanmean(F2MinusF3,1),nanstd(F2MinusF3,[],1)./sqrt(size(F2MinusF3,1)),'y',hax,hl);
elseif putTogetherThalData==1
    bigFres=bigFres(1,:);
    % Peak-norm. units
    for i=1:length(allUnits)
        allUnits(i).responses{1}=allUnits(i).responses{1}./max(max(allUnits(i).responses{1}));
    end
    % Take weighted average and variance across mice
    nmice=unique(fromMouse);
    nUnitsPerMouse=[14.8 7 7 5 9]*2;
    runningAv=zeros(size(allUnits(1).responses{1}));
    runningVar=zeros(size(allUnits(1).responses{1}));
    for i=1:length(nmice)
        scaleFac=nUnitsPerMouse(i)/sum(nUnitsPerMouse);
        subUnits=allUnits(fromMouse==nmice(i));
        rr=nan(size(subUnits(1).responses{1},1),size(subUnits(1).responses{1},2),length(subUnits));
        for j=1:length(subUnits)
            rr(:,:,j)=subUnits(j).responses{1};
        end
        runningAv=runningAv+scaleFac.*nanmean(rr,3);
        runningVar=runningVar+scaleFac.*nanvar(rr,[],3);
    end
    % Pick out F1,F2,F3
    fund=nan(1,length(freqs));
    F2=nan(1,length(freqs));
    F3=nan(1,length(freqs));
    fundVar=nan(1,length(freqs));
    F2Var=nan(1,length(freqs));
    F3Var=nan(1,length(freqs));
    for i=1:length(freqs)
        fund(i)=nanmean(runningAv(i,[find(bigFres>=freqs(i),1,'first') find(bigFres<=freqs(i),1,'last')]));
        F2(i)=nanmean(runningAv(i,[find(bigFres>=2*freqs(i),1,'first') find(bigFres<=2*freqs(i),1,'last')]));
        F3(i)=nanmean(runningAv(i,[find(bigFres>=3*freqs(i),1,'first') find(bigFres<=3*freqs(i),1,'last')]));
        fundVar(i)=nanmean(runningVar(i,[find(bigFres>=freqs(i),1,'first') find(bigFres<=freqs(i),1,'last')]));
        F2Var(i)=nanmean(runningVar(i,[find(bigFres>=2*freqs(i),1,'first') find(bigFres<=2*freqs(i),1,'last')]));
        F3Var(i)=nanmean(runningVar(i,[find(bigFres>=3*freqs(i),1,'first') find(bigFres<=3*freqs(i),1,'last')]));
    end
    % Plot results
    figure();
    hax=axes();
    hl=semilogx(freqs,fund,'Color','k');
    hold on;
    addErrBar(freqs,fund,sqrt(fundVar)./sqrt(sum(nUnitsPerMouse)),'y',hax,hl);
    hl=semilogx(freqs,F2-F3,'Color',[0.8 0.8 0.8]);
    hold on;
    addErrBar(freqs,F2-F3,sqrt(F2Var-F3Var)./sqrt(sum(nUnitsPerMouse)),'y',hax,hl);
    thalResults.fund=fund;
    thalResults.fundVar=fundVar;
else
    if doUnitsOneByOne==1
        for i=1:length(allUnits)
            [F1(i,:),F2(i,:),F3(i,:)]=plotFreqResponseWithHarms(allUnits(i).responses{1},bigFres);
            F2MinusF3(i,:)=F2(i,:)-F3(i,:);
        end
    else
        [F1,F2,F3]=plotFreqResponseWithHarms(allUnits,bigFres);
    end
    
    if integralAlign==1
        for i=1:size(F1,1)
            s=sum(F1(i,:)+F2(i,:)+F3(i,:));
            F1(i,:)=F1(i,:)./s;
            F2(i,:)=F2(i,:)./s;
            F3(i,:)=F3(i,:)./s;
            F2MinusF3(i,:)=F2MinusF3(i,:)./s;
        end
    end
end


    