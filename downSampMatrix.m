function newvec=downSampMatrix(matrix,n)

stepInds=1:n:size(matrix,2);
newvec=nan(size(matrix,1),length(stepInds));
% disp(length(stepInds));
for i=1:length(stepInds)
%     disp(i);
    if i==length(stepInds)
        newvec(:,i)=nanmean(matrix(:,stepInds(i):size(matrix,2)),2);
    else
        newvec(:,i)=nanmean(matrix(:,stepInds(i):stepInds(i+1)),2);
    end
end