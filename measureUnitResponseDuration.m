function allDurations=measureUnitResponseDuration(allUnitsTimeCourse,stimOn,spontWindow)

% Find duration from start of stim. to time point when response is back
% within noise (1 std dev) for at least thresh ms

dsf=100; % downSampFactor (in ms)
thresh=100; % time below noise mean + 1 std to call "end" of response
threshn=floor(thresh/dsf);

nodownx=allUnitsTimeCourse.x;
allDurations=nan(1,length(allUnitsTimeCourse.assigns));
for i=1:length(allUnitsTimeCourse.assigns)
    nodowny=allUnitsTimeCourse.y1(i,:);
    x=downSampAv(nodownx,dsf);
    y=downSampAv(nodowny,dsf);
    if sum(y==0)>0.5*length(y)
        x=downSampAv(nodownx,dsf*2);
        y=downSampAv(nodowny,dsf*2);
    end
    if sum(y==0)>0.5*length(y)
        x=downSampAv(nodownx,dsf*3);
        y=downSampAv(nodowny,dsf*3);
    end
    m=mean(y(x>=spontWindow(1) & x<spontWindow(2)));
    s=std(y(x>=spontWindow(1) & x<spontWindow(2)));
    isbelowx=x(find(y<m+s));
    isbelowx_ind=find(y<m+s);
    usebelowx=isbelowx(isbelowx>stimOn);
    usebelowx_ind=isbelowx_ind(isbelowx>stimOn);
    endOfResponse=nan;
    for j=1:length(usebelowx_ind)
        k=usebelowx_ind(j);
%         if ismember(k+1,usebelowx_ind) & ismember(k+2,usebelowx_ind) & ismember(k+3,usebelowx_ind) & ismember(k+4,usebelowx_ind) & ismember(k+5,usebelowx_ind) & ismember(k+6,usebelowx_ind)
        if ismember(k+1,usebelowx_ind) & ismember(k+2,usebelowx_ind) & ismember(k+3,usebelowx_ind)
            endOfResponse=usebelowx(j);
            break
        end
    end
    if isnan(endOfResponse)
        allDurations(i)=nan;
    else
        allDurations(i)=endOfResponse-stimOn;
    end
end