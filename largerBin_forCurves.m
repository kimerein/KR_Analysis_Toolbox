% for i=1:length(Ydata1)-1
%     newYdata1(i)=mean(Ydata1(i:i+1));
% end
    
a1=ypoints1;
a2=ypoints2;
x1=xpoints1;

n=30;
pickInds=1:n:length(a1)-(n-1);

b1=[];
b2=[];
c=[];

j=1;
for i=pickInds
    b1(j)=mean(a1(i:i+n-1));
    j=j+1;
end
j=1;
for i=pickInds
    b2(j)=mean(a2(i:i+n-1));
    j=j+1;
end
j=1;
for i=pickInds
    c(j)=mean(x1(i:i+n-1));
    j=j+1;
end
figure(); plot(c,b1,'Color','k'); hold on; plot(c,b2,'Color','r');

outYpoints1=b1;
outYpoints2=b2;
outXpoints1=c;
