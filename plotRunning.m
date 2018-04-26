function isRunning=plotRunning(encoder,LED,phys,i)

% % Code set for m113
% scaleFactor=0.9434;
% trialDuration=4;
% detectThresh=0.225;
% plotFigs=0;
% countThresh=3;

% % Code set for m116
% scaleFactor=0.9434;
% trialDuration=4;
% detectThresh=0.225;
% plotFigs=0;
% countThresh=2;

% % Code set for m206
% scaleFactor=0.9434;
% trialDuration=4;
% detectThresh=0.225;
% plotFigs=0;
% countThresh=2;

% Code set for m301
scaleFactor=0.7;
scaleFactor2=0.5;
trialDuration=10;
detectThresh=0.225;
plotFigs=1;
countThresh=2;
spontEnd=0.5;

times=linspace(0,trialDuration,size(LED,2));
y=encoder(i,:)-mean(encoder(1,linspace(0,trialDuration,size(LED,2))<spontEnd),2)-(scaleFactor)*(LED(i,:)-mean(LED(1,linspace(0,trialDuration,size(LED,2))<spontEnd),2))-(scaleFactor2)*(phys(i,:)-mean(phys(1,linspace(0,trialDuration,size(LED,2))<spontEnd),2));

if plotFigs==1
    figure();
    plot(times(5:end),encoder(i,5:end)-mean(encoder(1,linspace(0,trialDuration,size(LED,2))<spontEnd),2),'Color','k');
    hold on;
    plot(times(5:end),(scaleFactor)*(LED(i,5:end)-mean(LED(1,linspace(0,trialDuration,size(LED,2))<spontEnd),2))+(scaleFactor2)*(phys(i,5:end)-mean(phys(1,linspace(0,trialDuration,size(LED,2))<spontEnd),2)),'Color','b');
    
    figure(); 
    plot(times(5:end),y(5:end));
end
    
% bandpassed=bandPassLFP(y,1./(trialDuration/size(LED,2)),50,10000,0);
bandpassed=bandPassLFP(y,1./(trialDuration/size(LED,2)),100,200,0);

if plotFigs==1
    figure(); 
    plot(times(5:end),bandpassed(5:end));
end

fp=findpeaks(bandpassed(5:end));
if sum(fp>detectThresh)>=countThresh 
    isRunning=1;
else
    isRunning=0;
end