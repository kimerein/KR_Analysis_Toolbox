function [predictedCx,inputCx,outputCx,tToHalfMax,Pxy,f]=predictCorticalResponseGivenTau(tau_y_normInh,tau_y,xpoints,inputCx,outputCx)

spontWindow=[0 1];
normalize=0;
calcPhaseDelay=1;
% spontToSubtract=[0.6 1];
spontToSubtract=[1 1.05];
% peak=[1.2 1.3];
peak=[1.19 1.21];
beforeLED=[1.3 1.34];
afterLED=[1.46 1.5];
% top=[1.04+0.07 1.04+0.1];
% top=[1.2 1.3];
top=[1.19 1.21];
% startIncrease=1.02;
startIncrease=1.035;
stimOnset=1.04;
ledWindow=[0.6 2.9];

fillIn=1;
fills={[0.536 0.544]; [1.339 1.342]; [3.037 3.042]};
if fillIn==1
    for i=1:length(fills)
        currfill=fills{i};
        inputCx(xpoints>currfill(1) & xpoints<currfill(2))=ones(1,sum(xpoints>currfill(1) & xpoints<currfill(2)))*mean(inputCx(find(xpoints<currfill(1),10,'last')));
        outputCx(xpoints>currfill(1) & xpoints<currfill(2))=ones(1,sum(xpoints>currfill(1) & xpoints<currfill(2)))*mean(outputCx(find(xpoints<currfill(1),10,'last')));
    end
end

% First have to remove the tau under conditions of normal inhibition
% estimated_nsr=var(inputCx(xpoints>=spontWindow(1) & xpoints<=spontWindow(2)))/var(smooth(inputCx,100));
% estimated_nsr=100000000;
estimated_nsr=var(inputCx(xpoints>=spontWindow(1) & xpoints<=spontWindow(2)))/5;
% thal=deconvwnr(inputCx,tau_y_normInh,estimated_nsr);
thal=deconvwnr(inputCx,[zeros(1,length(tau_y_normInh)) tau_y_normInh],estimated_nsr);
% figure(); plot(thal);

predictedCx=conv(tau_y,thal);
predictedCx=predictedCx(1:length(inputCx));
scale=max(outputCx)./max(predictedCx);
predictedCx=predictedCx.*scale;

if normalize==1
    inputCx=inputCx-mean(inputCx(xpoints>=spontToSubtract(1) & xpoints<=spontToSubtract(2)));
    m=mean(inputCx(xpoints>=peak(1) & xpoints<=peak(2)));
    inputCx=inputCx./m;
    predictedCx=predictedCx-mean(predictedCx(xpoints>=spontToSubtract(1) & xpoints<=spontToSubtract(2)));
    m=mean(predictedCx(xpoints>=peak(1) & xpoints<=peak(2)));
    predictedCx=predictedCx./m;
    outputCx=outputCx-mean(outputCx(xpoints>=spontToSubtract(1) & xpoints<=spontToSubtract(2)));
    m=mean(outputCx(xpoints>=peak(1) & xpoints<=peak(2)));
    outputCx=outputCx./m;
%     predictedCx=predictedCx-mean(predictedCx(xpoints>=afterLED(1) & xpoints<=afterLED(2)));
%     predictedCx=predictedCx+mean(outputCx(xpoints>=afterLED(1) & xpoints<=afterLED(2)));
%     predictedCx=(mean(outputCx(xpoints>=beforeLED(1) & xpoints<=beforeLED(2)))/mean(predictedCx(xpoints>=beforeLED(1) & xpoints<=beforeLED(2)))).*predictedCx;
end

if calcPhaseDelay==1
    [Pxy,f]=cpsd(inputCx(xpoints>ledWindow(1) & xpoints<ledWindow(2)),outputCx(xpoints>ledWindow(1) & xpoints<ledWindow(2)),[],[],[],1000);
end
figure(); 
plot(xpoints(xpoints>=0.1),inputCx(xpoints>=0.1),'Color','k');
hold on; 
plot(xpoints(xpoints>=0.1),predictedCx(xpoints>=0.1),'Color','r');
plot(xpoints(xpoints>=0.1),outputCx(xpoints>=0.1),'Color','b');
[~,pt(1)]=max(inputCx(xpoints>=startIncrease & xpoints<=peak(2)));
useSub=xpoints(xpoints>=startIncrease & xpoints<=peak(2));
peakTimes(1)=useSub(pt(1)).*(useSub(2)-useSub(1));
[~,pt(2)]=max(predictedCx(xpoints>=startIncrease & xpoints<=peak(2)));
peakTimes(2)=useSub(pt(2)).*(useSub(2)-useSub(1));
[~,pt(3)]=max(outputCx(xpoints>=startIncrease & xpoints<=peak(2)));
peakTimes(3)=useSub(pt(3)).*(useSub(2)-useSub(1));
disp(['black ' 'red ' 'blue']);
disp(peakTimes);
figure();
plot(useSub,predictedCx(xpoints>=startIncrease & xpoints<=peak(2)));
hold on;
scatter(useSub(pt(2)),1);

topInput=mean(inputCx(xpoints>=top(1) & xpoints<=top(2)));
topPredicted=mean(predictedCx(xpoints>=top(1) & xpoints<=top(2)));
topOutput=mean(outputCx(xpoints>=top(1) & xpoints<=top(2)));
thresh{1}=topInput/2;
thresh{2}=topPredicted/2;
thresh{3}=topOutput/2;
fixedCurve{1}=inputCx(xpoints>=startIncrease & xpoints<=top(2));
fixedCurve{2}=predictedCx(xpoints>=startIncrease & xpoints<=top(2));
fixedCurve{3}=outputCx(xpoints>=startIncrease & xpoints<=top(2));
subX=xpoints(xpoints>=startIncrease & xpoints<=top(2));
tToHalfMax=zeros(1,3);
for i=1:length(fixedCurve)
    fixCurve=fixedCurve{i};
    th=thresh{i};
    fixCurve(fixCurve>th*2)=th*2;
    fixCurve(subX>=top(1))=th*2;
    downCross=find(fixCurve>th,1);
    upCross=find(fixCurve(end:-1:1)<th,1);
    upCross=length(fixCurve)-(upCross-1);
%     avCross=ceil(mean([downCross upCross]));
%     avCross=downCross;
%     crossTime=subX(avCross);
    crossTime=mean([subX(downCross) subX(upCross)]);
    avCross=find(subX>=crossTime,1);
    tToHalfMax(i)=crossTime-stimOnset;
    figure(); 
    plot(subX,fixCurve);
    hold on;
    scatter(crossTime,fixCurve(avCross),[],'r');
    line([top(1) top(2)],[th*2 th*2]);
end
disp(['black ' 'red ' 'blue']);
disp(tToHalfMax);