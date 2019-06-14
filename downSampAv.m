function newvec=downSampAv(vec,n)

if size(vec,1)>1
    newvec=zeros(length(1:n:length(vec)),1);
else 
    newvec=zeros(1,length(1:n:length(vec)));
end
stepInds=1:n:length(vec);
for i=1:length(stepInds)
    if i==length(stepInds)
        newvec(i)=nanmean(vec(stepInds(i):length(vec)));
    else
        newvec(i)=nanmean(vec(stepInds(i):stepInds(i+1)));
    end
end