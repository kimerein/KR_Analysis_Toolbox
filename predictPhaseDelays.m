function predictPhaseDelays(t,timeCourse1,timeCourse2)
% Make sure t is in seconds, not ms
% Could get timeCourse from toyModel

measurePhaseDelay=1;
transferFunction=[];

% result=getConvForFreq(10,1,t,timeCourse);

freq=1:1:60;
amps1=zeros(1,length(freq));
amps2=zeros(1,length(freq));
allPxy=zeros(1,length(freq));
timeDelays=zeros(1,length(freq));
for i=1:length(freq)
    [result1,x]=getConvForFreq(freq(i),0,t,timeCourse1,transferFunction);
    ma=max(result1);
    mi=min(result1);
    amps1(i)=ma-mi;
    [result2,x]=getConvForFreq(freq(i),0,t,timeCourse2,transferFunction);
    ma=max(result2);
    mi=min(result2);
    amps2(i)=ma-mi;
    if measurePhaseDelay==1
        [Pxy,fphase]=cpsd(result1,result2,[],[],[],1/(x(2)-x(1)));
        phase=angle(Pxy);
        allPxy(i)=mean(Pxy(fphase>freq(i)-0.1 & fphase<freq(i)+0.1));
        timeDelays(i)=mean(phase(fphase>freq(i)-0.1 & fphase<freq(i)+0.1))./(2*pi*freq(i));
    end
end
figure(); 
plot(freq,amps1,'Color','k');
hold on;
plot(freq,amps2,'Color','r');

if measurePhaseDelay==1
    figure(); 
    plot(freq,timeDelays);
end

end

function [out,x]=applyTransferFunction(x,in,transferFunction)
    th=transferFunction.thalamus;
    cx=transferFunction.cortex;
    out=zeros(size(in));
    for i=1:length(in)
        if in(i)<min(th)
            out(i)=0;
        elseif in(i)>max(th)
            out(i)=max(cx);
        else
            if isempty(find(th<in(i),1,'last'))
                out(i)=min(cx);
            elseif isempty(find(th>in(i),1,'first'))
                out(i)=max(cx);
            else
                out(i)=mean([cx(find(th<in(i),1,'last')) cx(find(th>in(i),1,'first'))]);
                if isnan(out(i))
                    disp('help');
                end
            end
        end
    end
end

function [output,x]=getConvForFreq(f,doPlot,t,timeCourse,transferFunction)
if isempty(timeCourse)
    x=0:0.0001:2; % in s
%     x=0:0.0001:2.0031; % in s
    net_tau=0.3; % in s
%     net_tau=0.02; % in s
%     net_tau=0.1; % in s
    shutoff=[ones(1,length(x(x>=0 & x<=0.003))) exp(-x./net_tau)];
%     shutoff=exp(-x./net_tau);
else
%     if max(t)<4
    if max(t)<20
        tsteps=t(2)-t(1);
%         x=[t max(t)+tsteps:tsteps:4];
        x=[t max(t)+tsteps:tsteps:20];
        timeCourse=[timeCourse zeros(1,length(max(t)+tsteps:tsteps:20))];
    else
        x=t; % in s
    end
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

ap=1;
input=ap*sin(2*pi*f.*x);
if ~isempty(transferFunction)
    [input,x]=applyTransferFunction(x,input,transferFunction);
end
output=conv(shutoff,input,'same');
if isempty(timeCourse)
    x=[x 2.0001:0.0001:2+0.003+0.0001];
end
end

