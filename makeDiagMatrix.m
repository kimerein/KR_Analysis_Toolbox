function makeDiagMatrix(diags)

n=8;
diagMatrix=zeros(n,15);
edges=linspace(0,1,n);
for i=1:size(diags,1)
    for j=1:size(diags,2)
        curry=diags(i,j);
        diagMatrix(find(edges>=curry,1,'first'),j)=diagMatrix(find(edges>=curry,1,'first'),j)+1;
    end
end

figure(); 
imagesc(diagMatrix);