function data=normalize_matrix_rows(xpoints,data,top,bottom)

% topInds=find(xpoints>3.25 & xpoints<3.3);
% bottomInds=find(xpoints>3.4 & xpoints<3.5);
topInds=find(xpoints>top(1) & xpoints<top(2));
bottomInds=find(xpoints>bottom(1) & xpoints<bottom(2));

% topInds=find(xpoints>1.19 & xpoints<1.2);
% bottomInds=find(xpoints>1.275 & xpoints<1.285);
% topInds=find(xpoints>1.0 & xpoints<1.2);
% bottomInds=find(xpoints>1.25 & xpoints<1.325);
% topInds=find(xpoints>1.022 & xpoints<1.03);
% bottomInds=find(xpoints>0.98 & xpoints<1);
% topInds=find(xpoints>1.01 & xpoints<1.04);
% bottomInds=find(xpoints>0.98 & xpoints<1);
% topInds=find(xpoints>1.365 & xpoints<1.375);
% bottomInds=find(xpoints>1.275 & xpoints<1.3);

% topInds=find(xpoints>1.35 & xpoints<1.38);
% bottomInds=find(xpoints>1.275 & xpoints<1.3);

sharedBaseline=mean(mean(data(:,bottomInds)));
for i=1:size(data,1)
    height=mean(data(i,topInds))-mean(data(i,bottomInds));
%     height=max(data(i,topInds))-mean(data(i,bottomInds));
    normFactor(i)=1/height;
end
for i=1:size(data,1)
    data(i,:)=(data(i,:)-mean(data(i,bottomInds)))*normFactor(i);
end