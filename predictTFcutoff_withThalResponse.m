function [giveOutput,freq]=predictTFcutoff_withThalResponse(t,timeCourse,thal_freq,thal_response,transferFunction)
% Make sure t is in seconds, not ms
% Could get timeCourse from toyModel

measurePhaseDelay=0;

% result=getConvForFreq(10,1,t,timeCourse);

freq=1:1:60;
amps=zeros(1,length(freq));
allPxy=zeros(1,length(freq));
timeDelays=zeros(1,length(freq));
for i=1:length(freq)
    tryInput=thal_response;
    tryInput(tryInput<0)=0;
    tryInput=sqrt(tryInput);
    [result,x]=getConvForFreq(freq(i),0,t,timeCourse,thal_freq,tryInput,transferFunction);
%     [result,x]=getConvForFreq(freq(i),0,t,timeCourse,thal_freq,sqrt(thal_response),transferFunction);
%     [result,x]=getConvForFreq(freq(i),0,t,timeCourse,thal_freq,thal_response,transferFunction);
    ma=max(result);
    mi=min(result);
%     [pks,loc]=findpeaks(result,'npeaks',1);
    amps(i)=ma-mi;
%     amps(i)=pks-min(result);
    if freq(i)==1
        figure(); 
        plot(x,result);
    end
    if measurePhaseDelay==1
        [Pxy,fphase]=cpsd(sin(2*pi*f.*x),result,[],[],[],1/(x(2)-x(1)));
        phase=angle(Pxy);
        allPxy(i)=mean(Pxy(fphase>freq(i)-0.1 & fphase<freq(i)+0.1));
        timeDelays(i)=mean(phase(fphase>freq(i)-0.1 & fphase<freq(i)+0.1))./(2*pi*fphase);
    end
end
figure(); 
plot(freq,amps.^2);
giveOutput=amps.^2;
% plot(freq,amps);
% giveOutput=amps;

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

function [output,x]=getConvForFreq(f,doPlot,t,timeCourse,thal_freq,thal_response,transferFunction)
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

% f=140; % in Hz
% Apply thalamic filter
ap=mean([thal_response(find(thal_freq<=f,1,'last')) thal_response(find(thal_freq>=f,1,'first'))]);
input=ap*sin(2*pi*f.*x);
if f==4
    figure();
    plot(x,input);
    title('Thalamic Frequency Response Gives Amplitude of This Sine');
end
if ~isempty(transferFunction)
    [input,x]=applyTransferFunction(x,input,transferFunction);
end
if f==4
    figure();
    plot(x,input);
    title('After Applying Transfer Function');
end
output=conv(shutoff,input,'same');
if isempty(timeCourse)
    x=[x 2.0001:0.0001:2+0.003+0.0001];
end
if f==4
    figure();
    plot(x,output);
    title('After Convolution with V1 Tau');
end

% figure(); 
% plot(output);
end

