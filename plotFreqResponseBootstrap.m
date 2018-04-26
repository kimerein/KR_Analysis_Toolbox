function [filledges,byCellResponses]=plotFreqResponseBootstrap(runningSU,bigFres,plotCell,c)

for i=1:length(runningSU(1).p)
    [fund,F2]=iterateThroughCells(runningSU,i,bigFres);
    if i==1
        byCellResponses=nan(size(fund,1),size(fund,2),length(runningSU(1).p));
        byCellResponses_F2=nan(size(fund,1),size(fund,2),length(runningSU(1).p));
        % each row is a cell, each column is a temp freq, each 3rd dim is
        % an iteration of the bootstrap
        byCellResponses(:,:,1)=fund;
        byCellResponses_F2(:,:,1)=F2;
    else
        byCellResponses(:,:,i)=fund;
        byCellResponses_F2(:,:,i)=F2;
    end
end

% c='b';
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% plotCell=1;
% figure(); 
semilogx(freqs,reshape(nanmean(byCellResponses(plotCell,:,:),3),1,length(freqs)),'Color',c,'LineWidth',2);
hold on;
% semilogx(freqs,reshape(prctile(byCellResponses(plotCell,:,:),2.5,3),1,length(freqs)),'Color',c,'LineWidth',1);
% semilogx(freqs,reshape(prctile(byCellResponses(plotCell,:,:),97.5,3),1,length(freqs)),'Color',c,'LineWidth',1);
% h=patch('XData',[freqs fliplr(freqs)],'YData',[reshape(prctile(byCellResponses(plotCell,:,:),2.5,3),1,length(freqs)) fliplr(reshape(prctile(byCellResponses(plotCell,:,:),97.5,3),1,length(freqs)))],'EdgeColor','none','FaceColor',c);
filledges=[reshape(prctile(byCellResponses(plotCell,:,:),2.5,3),1,length(freqs)) fliplr(reshape(prctile(byCellResponses(plotCell,:,:),97.5,3),1,length(freqs)))];
% Fix nans
try
    filledges=naninterp(filledges);
catch
end

h=patch('XData',[freqs fliplr(freqs)],'YData',filledges,'EdgeColor','none','FaceColor',c);
% h=fill([freqs fliplr(freqs)],filledges,c);
% set(h,'FaceColor',c);
% set(h,'facealpha',0.1);
% alpha(0.5);
semilogx(freqs,reshape(nanmean(byCellResponses(plotCell,:,:),3),1,length(freqs)),'Color',c,'LineWidth',2);

end

function [fund,F2]=iterateThroughCells(runningSU,currIt,bigFres)

for i=1:length(runningSU)
    [fund(i,:),F2(i,:)]=plotFreqResponseWithHarms_subFnctn(currIt,runningSU(i),bigFres,1,0,1);
end

end

function X = naninterp(X) 
% Interpolate over NaNs 
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'cubic'); 
return 
end

function [fund,F2,F3,r,rr,matrix,pp,inte]=plotFreqResponseWithHarms_subFnctn(currIt,r,bigFres,doAmp,plotFigs,subtractNonspec)

integralAlign=0;
peakAlign=0;

replaceAllZerosWithNan=1;

if doAmp==1
    if isfield(r,'responses')
        for i=1:length(r)
            r(i).responses{1}=sqrt(r(i).responses{currIt});
        end
    else
        r=sqrt(r);
    end
end

rr=[];
j=1;
if isfield(r,'responses')
    newR=zeros(size(r(1).responses{currIt}));
    for i=1:length(r)
        if integralAlign==1
            currR=r(i).responses{1}./sum(sum(r(i).responses{currIt}));
            currP=r(i).pStim{1}./sum(sum(r(i).responses{currIt}));
            inte(i)=sum(sum(r(i).responses{currIt}));
        elseif peakAlign==1
            currR=r(i).responses{1}./max(max(r(i).responses{currIt}));
            currP=r(i).pStim{1}./max(max(r(i).responses{currIt}));
        else
            currR=r(i).responses{currIt};
            currP=r(i).pStim{currIt};
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

end