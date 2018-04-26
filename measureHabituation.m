function [strong_resp,weak_resp,strong_base,weak_base]=measureHabituation(psth,useassigns)

removeRegression=1;

[~,~,~,~,~,tResp,tBase]=getCellOrdering(psth);

nTrials=size(psth.psths{1},1);
firstSet=floor(nTrials/2);

str=1; 
wea=1; 
for i=useassigns
    currb=tBase{i};
    currr=tResp{i};
    if nanmean(currr(1:firstSet))>1*nanmean(currb(1:firstSet))
        initStrongResp(str,:)=currr;
        initStrongResp_base(str,:)=currb;
        str=str+1;
    else
        initWeakResp(wea,:)=currr;
        initWeakResp_base(wea,:)=currb;
        wea=wea+1;
    end
    allResp(i,:)=currr;
    allResp_base(i,:)=currb;
end


if removeRegression==1
    figure(); plot(nanmean(initStrongResp,1)-nanmean(allResp,1),'Color','r'); hold on; plot(nanmean(initStrongResp_base,1)-nanmean(allResp_base,1),'Color','k');
    figure(); plot(nanmean(initWeakResp,1)-nanmean(allResp,1),'Color','r'); hold on; plot(nanmean(initWeakResp_base,1)-nanmean(allResp_base,1),'Color','k');
    
    strong_resp=nanmean(initStrongResp,1)-nanmean(allResp,1);
    weak_resp=nanmean(initWeakResp,1)-nanmean(allResp,1);
    strong_base=nanmean(initStrongResp_base,1)-nanmean(allResp_base,1);
    weak_base=nanmean(initWeakResp_base,1)-nanmean(allResp_base,1);
    
    figure(); plot(downSampAv(nanmean(initStrongResp,1)-nanmean(allResp,1),5),'Color','r'); hold on; plot(downSampAv(nanmean(initStrongResp_base,1)-nanmean(allResp_base,1),5),'Color','k');
    figure(); plot(downSampAv(nanmean(initWeakResp,1)-nanmean(allResp,1),5),'Color','r'); hold on; plot(downSampAv(nanmean(initWeakResp_base,1)-nanmean(allResp_base,1),5),'Color','k');
else
    figure(); plot(nanmean(initStrongResp,1),'Color','r'); hold on; plot(nanmean(initStrongResp_base,1),'Color','k');
    figure(); plot(nanmean(initWeakResp,1),'Color','r'); hold on; plot(nanmean(initWeakResp_base,1),'Color','k');
    
    strong_resp=nanmean(initStrongResp,1);
    weak_resp=nanmean(initWeakResp,1);
    strong_base=nanmean(initStrongResp_base,1);
    weak_base=nanmean(initWeakResp_base,1);
    
    figure(); plot(downSampAv(nanmean(initStrongResp,1),5),'Color','r'); hold on; plot(downSampAv(nanmean(initStrongResp_base,1),5),'Color','k');
    figure(); plot(downSampAv(nanmean(initWeakResp,1),5),'Color','r'); hold on; plot(downSampAv(nanmean(initWeakResp_base,1),5),'Color','k');
end