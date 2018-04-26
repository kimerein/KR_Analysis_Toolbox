function predictTFcutoff(t,timeCourse)
% Make sure t is in seconds, not ms
% Could get timeCourse from toyModel

% result=getConvForFreq(10,1,t,timeCourse);

freq=1:1:100;
amps=zeros(1,length(freq));
for i=1:length(freq)
    [result,x]=getConvForFreq(freq(i),0,t,timeCourse);
    ma=max(result);
    mi=min(result);
    amps(i)=ma-mi;
    if freq(i)==10
        figure(); 
        plot(x,result);
    end
end
figure(); 
plot(freq,amps);
end

function [output,x]=getConvForFreq(f,doPlot,t,timeCourse)
if isempty(timeCourse)
    x=0:0.0001:2; % in s
    net_tau=0.012; % in s
    shutoff=[ones(1,length(x(x>=0 & x<=0.003))) exp(-x./net_tau)];
else
    x=t; % in s
    net_tau=0.012; % in s
    % Scale starting value to 1
    ma=max(timeCourse);
    scaleFactor=1/ma;
    shutoff=timeCourse.*scaleFactor;
end

if doPlot==1
    figure(); 
    if isempty(timeCourse)
        plot([x 2.0001:0.0001:2+0.003+0.0001],shutoff);
    else
        plot(x,shutoff);
    end
end

% f=140; % in Hz
input=sin(2*pi*f.*x);
output=conv(shutoff,input,'same');
if isempty(timeCourse)
    x=[x 2.0001:0.0001:2+0.003+0.0001];
end

% figure(); 
% plot(output);
end

