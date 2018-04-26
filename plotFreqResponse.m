function [diags,allY,allY_low,allY_high,insteadMatrix]=plotFreqResponse(diags,insteadMatrix,grps,perc,perc_low,perc_high,showFigs,bigFres,peakNorm,freqs)

% peakNorm=1;
integralAlign=1;
doCutoff=0;
percCutoff=40;
useLPindex=0;
% freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60]; % for Chirps

if ~isempty(insteadMatrix)
    if size(insteadMatrix,2)==15
        for i=1:15
            for j=1:size(insteadMatrix,3)
                diags(j,i)=insteadMatrix(i,i,j);
            end
        end
    else
        for i=1:length(freqs)
            for j=1:size(insteadMatrix,3)
                diags(j,i)=nanmean(insteadMatrix(i,[find(bigFres(1,:)>=freqs(i),1,'first') find(bigFres(1,:)<=freqs(i),1,'last')],j),2);
            end
        end
    end
    if peakNorm==1
        for j=1:size(insteadMatrix,3)
            insteadMatrix(:,:,j)=insteadMatrix(:,:,j)./max(diags(j,:));
        end
    end
end

if peakNorm==1
    for i=1:size(diags,1)
        diags(i,:)=diags(i,:)./max(diags(i,:));
    end
end


if integralAlign==1
    runningSum=0;
    for i=1:size(diags,1)
        diags(i,:)=diags(i,:)./sum(diags(i,:));
        runningSum=runningSum+sum(diags(i,:));
    end
    diags=diags*runningSum*size(diags,1);
end

% figure();
% semilogx(freqs,diags);

if useLPindex==1
    for i=1:size(diags,1)
        LPindex(i)=nanmean(diags(i,1:8))./nanmean(diags(i,9:15));
    end
    if showFigs==1
        figure(); 
        hax=axes();
    end
    Y=prctile(LPindex',perc,1);
    if showFigs==1
        hl=semilogx(freqs,nanmean(diags(LPindex>=Y,:),1));
        hold on;
    end
    Y_low=prctile(LPindex',perc_low,1);
    Y_high=prctile(LPindex',perc_high,1);
    Y_low_temp=nanmean(diags(LPindex>=Y_low,:),1);
    Y_high_temp=nanmean(diags(LPindex>=Y_high,:),1);
    if showFigs==1
        semilogx(freqs,Y_low_temp);
        semilogx(freqs,Y_high_temp);
    end
    allY=nanmean(diags(LPindex>=Y,:));
    allY_low=Y_low_temp;
    allY_high=Y_high_temp;
    return
end

if doCutoff==1
    for i=1:size(diags,2)
        temp=sort(diags(:,i));
%         tempdiags(1:length(temp(ceil((percCutoff/100)*length(temp)):end)),i)=temp(ceil((percCutoff/100)*length(temp)):end);
        tempdiags(1:length(temp(1:ceil((percCutoff/100)*length(temp)))),i)=temp(1:ceil((percCutoff/100)*length(temp)));
    end
    diags=tempdiags;
end

if ~isempty(grps)
    g=unique(grps);
    if showFigs==1
        figure();
        hax=axes();
    end
    for i=1:length(g)
        Y=prctile(diags(grps==g(i),:),perc,1);
        if showFigs==1
            hl=semilogx(freqs,Y);
            hold on;
        end
        Y_low=prctile(diags(grps==g(i),:),perc_low,1);
        Y_high=prctile(diags(grps==g(i),:),perc_high,1);
        Y_low=Y-Y_low;
        Y_high=Y_high-Y;
        if showFigs==1
            addErrBar_asymmetric(freqs,Y,Y_low,Y_high,'y',hax,hl);
            hold all;
        end
        allY(i,:)=Y;
        allY_low(i,:)=Y_low;
        allY_high(i,:)=Y_high;
    end
else
    if showFigs==1
        figure(); 
        hax=axes();
    end
    Y=prctile(diags,perc,1);
    if showFigs==1
        hl=semilogx(freqs,Y);
        hold on;
    end
    Y_low=prctile(diags,perc_low,1);
    Y_high=prctile(diags,perc_high,1);
    Y_low=Y-Y_low;
    Y_high=Y_high-Y;
    if showFigs==1
        addErrBar_asymmetric(freqs,Y,Y_low,Y_high,'y',hax,hl);
    end
    allY=Y;
    allY_low=Y_low;
    allY_high=Y_high;
end
end

% function diags=plotFreqResponse(diags,insteadMatrix,grps)
% 
% peakNorm=0;
% integralAlign=0;
% freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% 
% if ~isempty(insteadMatrix)
%     for i=1:15
%         for j=1:size(insteadMatrix,3)
%             diags(j,i)=insteadMatrix(i,i,j);
%         end
%     end
% end
% 
% if peakNorm==1
%     for i=1:size(diags,1)
%         diags(i,:)=diags(i,:)./max(diags(i,:));
% %         diags(i,:)=diags(i,:)./(diags(i,5));
%     end
% end
% 
% if integralAlign==1
%     runningSum=0;
%     for i=1:size(diags,1)
%         diags(i,:)=diags(i,:)./sum(diags(i,:));
%         runningSum=runningSum+sum(diags(i,:));
%     end
%     diags=diags*runningSum*size(diags,1);
% end
% 
% % figure();
% % semilogx(freqs,diags);
% 
% if ~isempty(grps)
%     g=unique(grps);
%     figure();
%     hax=axes();
%     for i=1:length(g)
%         hl=semilogx(freqs,nanmean(diags(grps==g(i),:),1));
%         hold on;
%         addErrBar(freqs,nanmean(diags(grps==g(i),:),1),nanstd(diags(grps==g(i),:),[],1)./sqrt(size(diags(grps==g(i),:),1)),'y',hax,hl);
%         hold all;
%     end
% else
%     figure(); 
%     hax=axes();
%     hl=semilogx(freqs,nanmean(diags,1));
%     hold on;
%     addErrBar(freqs,nanmean(diags,1),nanstd(diags,[],1)./sqrt(size(diags,1)),'y',hax,hl);
% end
