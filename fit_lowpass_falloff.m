function p=fit_lowpass_falloff(diags_thal,diags_cx,freqs,fitFreqs,perc,percLow,percHigh)

[~,thalCurve]=plotFreqResponse_nestedFunction(diags_thal,[],[],perc,percLow,percHigh,0,[]);
[~,cxCurve]=plotFreqResponse_nestedFunction(diags_cx,[],[],perc,percLow,percHigh,0,[]);

thalCurveZero=thalCurve(fitFreqs);
cxCurveZero=cxCurve(fitFreqs);
log_y=log(cxCurveZero./thalCurveZero);
log_x=freqs(fitFreqs);

figure(); plot(log(freqs),log(cxCurve./thalCurve));

p=polyfit(log_x',log_y',1);

figure();
plot(log_x,log_y);
hold on;
plot(log_x,log_x.*p(1)+p(2),'Color','r');

end

function [diags,allY,allY_low,allY_high,insteadMatrix]=plotFreqResponse_nestedFunction(diags,insteadMatrix,grps,perc,perc_low,perc_high,showFigs,bigFres)

peakNorm=1;
integralAlign=0;
doCutoff=0;
percCutoff=60;
useLPindex=0;
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

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
