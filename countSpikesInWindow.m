function [sumInWindow,ledTrials,xout,nout,noutOff]=countSpikesInWindow(x,n,timeWindow,ledOn,ledOff,ledTrials)

xout=[];
nout=[];
noutOff=[];

sumInWindow=zeros(1,size(n,1));
for i=1:size(n,1)
    sumInWindow(i)=sum(n(i,x>=timeWindow(1) & x<=timeWindow(2)));
end

% figure(); 
% [nout,xout]=hist(sumInWindow(ledTrials==ledOn),10);
% plot(xout,nout,'Color','r');
% hold on; 
% [noutOff,xout]=hist(sumInWindow(ledTrials==ledOff),10);
% plot(xout,noutOff,'Color','k');