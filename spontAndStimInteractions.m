function spontAndStimInteractions(spikes)

% Stim. constant, spont. changes
windowWidth=0.05;
windowOffset=0.025;
windowTimewindow=[0 1];
ledVals=[0];
bins=50;
divideData=3;

a=windowTimewindow(1);
b=windowTimewindow(2);
i=1;
while a<b
    if a+windowWidth>b
        break
    end
    window{i}=[a a+windowWidth];
    led{i}=ledVals;
    a=a+windowOffset;
    i=i+1;
end

alln=[];
for i=1:length(led)
    if isnan(ledVals)
        [m,s,n]=calcMeanAndStdDuringWindow(spikes,window{i});
    else
        [m,s,n]=calcMeanAndStdDuringWindow(filtspikes(spikes,0,'led',led{i}),window{i});
    end
    batchSize=floor(length(n)/3);
    n1=n(1:batchSize);
    n2=n(batchSize+1:batchSize*2);
    n3=n(batchSize*2+1:end);
    alln=[alln; mean(n1); mean(n2); mean(n3)];
    disp(i);
end

[heights,centers]=hist(alln,bins);
centers=[centers(1)-(centers(2)-centers(1)) centers centers(end)+(centers(2)-centers(1))];
heights=[0 heights 0];
figure(); 
hist(alln,bins);
hold on;
plot(centers,heights);
fitted=fit(centers',heights','gauss2');
% fitx=centers(1)-(centers(2)-centers(1)):0.1:centers(end)+(centers(2)-centers(1));
fitx=0:0.1:centers(end)+(centers(2)-centers(1));
y1=fitted.a1.*exp(-((fitx-fitted.b1)./fitted.c1).^2);
plot(fitx,y1,'Color','r');
y2=fitted.a2.*exp(-((fitx-fitted.b2)./fitted.c2).^2);
plot(fitx,y2,'Color','g');
y3=y1+y2;
plot(fitx,y3,'Color','y');
interse=find(y2>y1,1);
disp('x value at intersection');
disp(fitx(interse));

% Spont. constant, stim. changes
windowTimewindow=[0 1];
ledVals=[0];
bins=50;

alln=[];
if isnan(ledVals)
    [m,s,n]=calcMeanAndStdDuringWindow(spikes,windowTimewindow);
else
    [m,s,n]=calcMeanAndStdDuringWindow(filtspikes(spikes,0,'led',ledVals),windowTimewindow);
end
alln=n;

[heights,centers]=hist(alln,bins);
centers=[centers(1)-(centers(2)-centers(1)) centers centers(end)+(centers(2)-centers(1))];
heights=[0 heights 0];
figure(); 
hist(alln,bins);
hold on;
plot(centers,heights);
fitted=fit(centers',heights','gauss2');
% fitx=centers(1)-(centers(2)-centers(1)):0.1:centers(end)+(centers(2)-centers(1));
fitx=0:0.1:centers(end)+(centers(2)-centers(1));
y1=fitted.a1.*exp(-((fitx-fitted.b1)./fitted.c1).^2);
plot(fitx,y1,'Color','r');
y2=fitted.a2.*exp(-((fitx-fitted.b2)./fitted.c2).^2);
plot(fitx,y2,'Color','g');
y3=y1+y2;
plot(fitx,y3,'Color','y');
interse=find(y2>y1,1);
disp('x value at intersection');
disp(fitx(interse));


