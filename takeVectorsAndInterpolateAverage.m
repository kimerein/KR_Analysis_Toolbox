function [outx,outy,outy_std,ally]=takeVectorsAndInterpolateAverage(x,y,usex)

% x and y are cell arrays with corresponding elements of matching-length
% vectors giving x and y data sampled at various time points
% usex is a vector indicating which sample time points to use for interp

newy=zeros(length(x),length(usex));
for i=1:length(x)
    currx=x{i};
    curry=y{i};
    newys=interp1(currx,curry,usex,'linear',nan);
    newy(i,:)=newys;
end

outx=usex;
outy=nanmean(newy,1);
outy_std=nanstd(newy,0,1);
ally=newy;
