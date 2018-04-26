function [h_matrix,p_matrix]=makeSignificanceMatrix(data,alpha,tail)

h_matrix=zeros(size(data,2),size(data,2));
p_matrix=zeros(size(data,2),size(data,2));

for j=1:15
    for i=j+1:size(data,2)
%         disp(i);
%         disp(j);
        [h,p]=ttest(data(:,j),data(:,i),alpha,tail);
        h_matrix(j,i)=h;
        p_matrix(j,i)=p;
    end
end