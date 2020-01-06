function out=concatPSTHsFromSameExpt(curr,next)

f=fieldnames(curr);
for i=1:length(f)
    if strcmp(f{i},'t')
        continue
    end
    temp=cell(length(curr.(f{i}))+length(next.(f{i})),1);
    temp(1:length(curr.(f{i})))=curr.(f{i});
    temp(length(curr.(f{i}))+1:end)=next.(f{i});
    out.(f{i})=temp;
end
out.t=curr.t;

