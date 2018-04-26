
Es=zeros(length(0.005:0.01:1),1);
Ei=zeros(length(0.005:0.01:1),1);
ps=0.005:0.01:1;
for j=1:length(ps)
p_ee=ps(j);
p_ei=1;
p_inh=0.1;
p_exc=1-p_inh;
n_tot=2000;

% runningSum=0;
% for i=1:10000
%     runningSum=runningSum+i*((1-p_ee)^(i-1))*p_ee*p_exc^i;
% end
% disp(runningSum);

% runningSum=0;
% for i=1:10000
%     runningSum=runningSum+i*((1-p_ee)^(i-1))*p_ei*p_exc^(i-1)*p_inh;
% end
% disp(runningSum);

a=n_tot*p_ee*(1-p_ee);
n=2.433;
n=2000;
E_distance=(1/(1-a^n))*((1-n*a^n-a^n+a*n*a^n)/(1-a));
% disp(E_distance);
E_distance=(1-(1-p_ee)^(n+1))/p_ee-(n+1)*(1-p_ee)^n;
Es(j)=E_distance;

num=0;
den=0;
n=2.433;
for k=1:n
    num=num+(1-p_ee)^(k-1)*(n_tot*p_ee)^k*p_ee*k;
    den=den+(1-p_ee)^(k-1)*(n_tot*p_ee)^k*p_ee;
end
Es(j)=E_distance;

num=0;
den=0;
n=2.08;
n=100;
n_i=70;
for k=1:n
    num=num+(1-p_ee)^(k-1)*(n_tot*p_ee)^(k-1)*n_i*k;
    den=den+(1-p_ee)^(k-1)*(n_tot*p_ee)^(k-1)*n_i;
end
Ei(j)=num/den;


end

figure();
plot(ps,Es);
hold on;
plot(ps,Ei,'Color','r');