function [outx,result]=predictCxResponse_fromThalResponse(t,timeCourse,th_response,transferFunction)
% Make sure t is in seconds, not ms
% Could get timeCourse from toyModel

[result,x]=getConvForFreq(th_response,0,t,timeCourse,transferFunction);
figure(); 
outx=0:x(2)-x(1):(x(2)-x(1))*(length(result)-1);
plot(outx,result);
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

function [output,x]=getConvForFreq(input,doPlot,t,timeCourse,transferFunction)
if isempty(timeCourse)
    x=0:0.0001:2; % in s
%     x=0:0.0001:2.0031; % in s
%     net_tau=0.012; % in s
    net_tau=0.02; % in s
    shutoff=[ones(1,length(x(x>=0 & x<=0.003))) exp(-x./net_tau)];
%     shutoff=exp(-x./net_tau);
else
    if max(t)<4
        tsteps=t(2)-t(1);
        x=[t max(t)+tsteps:tsteps:4];
        timeCourse=[timeCourse zeros(1,length(max(t)+tsteps:tsteps:4))];
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

% Apply thalamic filter
if doPlot==1
    figure();
    plot(x,input);
    title('Thalamic Response');
end
    
if ~isempty(transferFunction)
    [input,x]=applyTransferFunction(x,input,transferFunction);
end

if doPlot==1
    figure();
    plot(x,input);
    title('After Applying Transfer Function');
end
    
% output=conv(shutoff,input,'same');
output=conv(shutoff,input);
if isempty(timeCourse)
    x=[x 2.0001:0.0001:2+0.003+0.0001];
end

% figure();
% plot(0:x(2)-x(1):(x(2)-x(1))*(length(output)-1),output);
% title('After Convolution with V1 Tau');

end

