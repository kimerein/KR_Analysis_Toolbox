function [m2A,m2B]=solveForFrequencyDoubling2(V_F1,V_F2,d_F1,d_F2,F_CDF)

tryScaleFacs=0.001:0.001:2;
err1=zeros(1,length(tryScaleFacs));
err2=zeros(1,length(tryScaleFacs));
for i=1:length(tryScaleFacs)
    [err1(i),err2(i)]=calcSubfunction(V_F1,V_F2,d_F1,d_F2,F_CDF,tryScaleFacs(i),0);
end

figure(); 
plot(tryScaleFacs,err1,'Color','k');
hold on; 
plot(tryScaleFacs,err2,'Color','r');

[~,~,m2A,m2B]=calcSubfunction(V_F1,V_F2,d_F1,d_F2,F_CDF,1.007,1);

end

function [err1,err2,model2A,model2B]=calcSubfunction(V_F1,V_F2,d_F1,d_F2,F_CDF,scaleFac,showFigs)
% scaleFac=0.85;
notBelowZero=1;
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

c_0=mean(V_F1)./(scaleFac.*mean(d_F1));
% Model 2A - no filtering by CDF
D_fromeq1=1-(V_F1./(c_0.*d_F1)); % D as a function of f
D_fromeq2=(V_F2-c_0.*d_F2)./(c_0.*d_F1); % D as a function of f
if notBelowZero==1
    D_fromeq1(D_fromeq1<0)=0;
    D_fromeq2(D_fromeq2<0)=0;
end

c_0=mean(V_F1)./(scaleFac.*mean(d_F1.*F_CDF));
% Model 2B - filtering by CDF
D_fromeq1_2B=1-(V_F1./(c_0.*F_CDF.*d_F1)); % D as a function of f
D_fromeq2_2B=(V_F2-c_0.*F_CDF.*d_F2)./(c_0.*F_CDF.*d_F1); % D as a function of f
if notBelowZero==1
    D_fromeq1_2B(D_fromeq1_2B<0)=0;
    D_fromeq2_2B(D_fromeq2_2B<0)=0;
end

if showFigs==1
    figure();
    semilogx(freqs,D_fromeq1,'-k');
    hold on;
    semilogx(freqs,D_fromeq2,'-r');
    title('Model 2A no LP filtering');
    axis([1 60 0 1]);
end
model2A.eq1=D_fromeq1;
model2A.eq2=D_fromeq2;
err1=sum(abs(D_fromeq1(1:12)-D_fromeq2(1:12)));
% disp(err1);

if showFigs==1
    figure();
    semilogx(freqs,D_fromeq1_2B,'-k');
    hold on;
    semilogx(freqs,D_fromeq2_2B,'-r');
    title('Model 2B w LP filtering');
    axis([1 60 0 1]);
end
model2B.eq1=D_fromeq1_2B;
model2B.eq2=D_fromeq2_2B;
err2=sum(abs(D_fromeq1_2B(1:12)-D_fromeq2_2B(1:12)));
% disp(err2);
end

