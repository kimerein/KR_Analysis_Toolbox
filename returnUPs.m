function UPs=returnUPs(xpoints,ypoints,UP_thresh)

% Return UPs
in_UP=0;
UP_starts=[];
UP_ends=[];
for i=1:length(xpoints)
    if in_UP==0
        if ypoints(i)>=UP_thresh
            in_UP=1;
            UP_starts=[UP_starts; xpoints(i)];
        end
    else
        if ypoints(i)<UP_thresh
            in_UP=0;
%             UP_ends=[UP_ends; xpoints(i)];
            UP_ends=[UP_ends; xpoints(i-1)];
        end
    end
end
if in_UP==1
    UP_ends=[UP_ends; xpoints(end)];
end
UPs=[UP_starts UP_ends];