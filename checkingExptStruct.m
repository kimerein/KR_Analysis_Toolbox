% if sum(expt.files.triggers)~=length(expt.sweeps.fileInd)
%     disp('Trials mismatched. Need to fix.');
% else
%     disp('Trials matched in number.');
% end
% 
% for i=1:length(expt.files.triggers)
%     
%     
b=[];
j=1;
for i=1:72
    s=0;
    while a(j)==i
        s=s+1;
        j=j+1;
        if j==length(a)
            break
        end
    end
    disp(s);
    b=[b s];
end