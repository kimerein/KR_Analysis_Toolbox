function [fund,F2,F3,r,rr,matrix,pp,inte]=plotFreqResponseWithHarms(r,bigFres,doAmp,plotFigs,subtractNonspec)

integralAlign=1;
peakAlign=0;
inte=[];

replaceAllZerosWithNan=1;

if doAmp==1
    if isfield(r,'responses')
        for i=1:length(r)
            r(i).responses{1}=r(i).responses{1}/(10^4);
            r(i).responses{1}=sqrt(r(i).responses{1}*100);
        end
    else
        r=r/(10^4);
        r=sqrt(r*100);
    end
else
    if isfield(r,'responses')
        for i=1:length(r)
            r(i).responses{1}=r(i).responses{1}/(10^4);
        end
    else
        r=r/(10^4);
    end
end

rr=[];
j=1;
if isfield(r,'responses')
    newR=zeros(size(r(1).responses{1}));
    for i=1:length(r)
        if integralAlign==1
            currR=r(i).responses{1}./sum(sum(r(i).responses{1}));
            currP=r(i).pStim{1}./sum(sum(r(i).responses{1}));
            inte(i)=sum(sum(r(i).responses{1}));
        elseif peakAlign==1
            currR=r(i).responses{1}./max(max(r(i).responses{1}));
            currP=r(i).pStim{1}./max(max(r(i).responses{1}));
        else
            currR=r(i).responses{1};
            currP=r(i).pStim{1};
        end
%         if max(currR)>100
%             newR=newR+currR;
%             rr(:,:,j)=currR;
%             j=j+1;
%         end 
        newR=newR+currR;
        rr(:,:,j)=currR;
        if size(currP,1)~=15
            currP=currP(1:15,1:15);
        end
        pp(:,:,j)=currP;
        j=j+1;
    end
%     disp(j);
    r=newR;
end

matrix=r(:,bigFres(1,:)>=0.5 & bigFres(1,:)<=60.5);
bigFres2=bigFres(1,bigFres(1,:)>=0.5 & bigFres(1,:)<=60.5);
bigFres=bigFres(1,:);

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

if replaceAllZerosWithNan==1
    matrix(matrix==0)=nan;
    r(r==0)=nan;
    rr(rr==0)=nan;
end

if subtractNonspec==1
    for i=1:size(matrix,1)
        currfre=freqs(i);
%     for i=1:size(bigFres2,1)
%         currfre=bigFres2(i);
%         matrix(:,i)=matrix(:,i)-nanmean(matrix(find(freqs>currfre,1,'first')+1:end,i));
        matrix(i,:)=matrix(i,:)-nanmean(r(i,find(bigFres>currfre*2,1,'first')+1:end));
        r(i,:)=r(i,:)-nanmean(r(i,find(bigFres>currfre*2,1,'first')+1:end));
    end
end

if plotFigs==1
    figure();
    imagesc(r);
end
    
fund=nan(1,length(freqs));
F2=nan(1,length(freqs));
F3=nan(1,length(freqs));
for i=1:length(freqs)
    fund(i)=nanmean(matrix(i,[find(bigFres2>=freqs(i),1,'first') find(bigFres2<=freqs(i),1,'last')]));
    F2(i)=nanmean(matrix(i,[find(bigFres2>=2*freqs(i),1,'first') find(bigFres2<=2*freqs(i),1,'last')]));
    F3(i)=nanmean(matrix(i,[find(bigFres2>=3*freqs(i),1,'first') find(bigFres2<=3*freqs(i),1,'last')]));
end

frr=nan(size(rr,3),length(freqs));
if ~isempty(rr)
    for j=1:size(rr,3)
        for i=1:length(freqs)
            frr(j,i)=nanmean(rr(i,[find(bigFres>=freqs(i),1,'first') find(bigFres<=freqs(i),1,'last')],j));
        end
    end
end

if plotFigs==1
%     figure();
%     semilogx(freqs,frr);
    
    figure();
    semilogx(freqs,fund,'Color','k');
    hold on;
    semilogx(freqs,F2,'Color','r');
%     semilogx(freqs,F3,'Color','g');
%     semilogx(freqs,fund+F2,'Color','b');
    
%     figure();
%     semilogx(freqs,fund,'Color','k');
%     hold on;
%     F2minusNoise=F2-F3;
%     F2minusNoise(F2minusNoise<0)=0;
%     semilogx(freqs,F2minusNoise,'Color',[0.8 0.8 0.8]);
%     semilogx(freqs,fund+F2minusNoise,'Color','b');
end