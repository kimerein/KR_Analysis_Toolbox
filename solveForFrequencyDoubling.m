function solveForFrequencyDoubling(V_F1,V_F2,d_F1,d_F2,F_CDF)

% c_0_test=0.01:0.05:1000;

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

% devs_model2A=nan(size(c_0_test));
% devs_model2B=nan(size(c_0_test));
% for i=1:length(c_0_test)
%     c_0=c_0_test(i);
%     devs_model2A(i)=sum(V_F1-c_0.*1.*d_F1);
%     devs_model2B(i)=sum(V_F1-c_0.*1.*F_CDF.*d_F1);
% end
  
c_0_lessthan_F1=V_F1./d_F1;
c_0_lessthan_F2=V_F2./d_F2;
% c_0_test_end=min([c_0_lessthan_F1 c_0_lessthan_F2]);
% devs_model2A(c_0_test>min(c_0_lessthan_F1))=nan;
% c_0_greaterthan_F1=(1+V_F1)./d_F1;
c_0_greaterthan_F2=V_F2./(d_F2+d_F1);
% devs_model2A(c_0_test<max(c_0_greaterthan_F2))=nan;
% c_0_test_start=max(c_0_greaterthan_F2);

c_0_test_start=max([c_0_lessthan_F1 c_0_lessthan_F2]);
c_0_test_end=min(c_0_greaterthan_F2);

c_0_test=c_0_test_start:(c_0_test_end-c_0_test_start)/1000:c_0_test_end;
tryD=0:0.001:1;
atEachFreq=cell(1,length(freqs)-3);
for k=1:length(freqs)-3
    disp(k);
    devs_model2A=nan(length(c_0_test),length(tryD));
    devs_model2A_F2=nan(length(c_0_test),length(tryD));
    for i=1:length(c_0_test)
        for j=1:length(tryD)
            c_0=c_0_test(i);
            devs_model2A(i,j)=(V_F1(k)-c_0.*(1-tryD(j)).*d_F1(k));
            devs_model2A_F2(i,j)=(V_F2(k)-c_0.*d_F2(k)-c_0.*tryD(j).*d_F1(k));
        end
    end
%     atEachFreq{k}=min(devs_model2A,[],2)+min(devs_model2A_F2,[],2);
    atEachFreq{k}=min(devs_model2A,[],2);
end
allErrs=zeros(size(atEachFreq{1}));
for k=1:length(freqs)-3
    allErrs=allErrs+atEachFreq{k};
end

[mi,mind]=nanmin(allErrs);
c_0=c_0_test(mind);
% c_0=mean([c_0_lessthan_F1 c_0_lessthan_F2 c_0_greaterthan_F2]);
c_0=mean(V_F1)./(1.*mean(d_F1));
% Model 2A - no filtering by CDF
D_fromeq1=1-(V_F1./(c_0.*d_F1)); % D as a function of f
D_fromeq2=(V_F2-c_0.*d_F2)./(c_0.*d_F1); % D as a function of f
% V_F1_pred=c_0.*(1-D_fromeq2).*d_F1;

c_0_lessthan_F1=V_F1./(d_F1.*F_CDF);
c_0_lessthan_F2=V_F2./(d_F2.*F_CDF);
% c_0_test_end=min([c_0_lessthan_F1 c_0_lessthan_F2]);

% devs_model2B(c_0_test>min(c_0_lessthan_F1))=nan;
% c_0_greaterthan_F1=(1+V_F1)./(d_F1.*F_CDF);
c_0_greaterthan_F2=V_F2./((d_F2+d_F1).*F_CDF);
% devs_model2B(c_0_test<max(c_0_greaterthan_F2))=nan;
% c_0_test_start=max([c_0_lessthan_F1 c_0_lessthan_F2 c_0_greaterthan_F2]);
% c_0_test_start=max(c_0_greaterthan_F2);

c_0_test_start=max([c_0_lessthan_F1 c_0_lessthan_F2]);
c_0_test_end=min(c_0_greaterthan_F2);

% c_0_test=c_0_test_start:(c_0_test_end-c_0_test_start)/1000:c_0_test_end;
c_0_test=0:1:10000;
canUse=zeros(size(c_0_test));
% c_0_test=c_0_test_start:0.1:1000;
tryD=0:0.001:1;
atEachFreq=cell(1,length(freqs)-3);
for i=1:length(c_0_test)
    c_0=c_0_test(i);
    D_fromeq1_2B=1-(V_F1./(c_0.*F_CDF.*d_F1)); 
    D_fromeq2_2B=(V_F2-c_0.*F_CDF.*d_F2)./(c_0.*F_CDF.*d_F1); 
    canUse(i)=all(D_fromeq1_2B>=0) & all(D_fromeq2_2B>=0) & all(D_fromeq1_2B<=1) & all(D_fromeq2_2B<=1);
end
for k=1:length(freqs)-3
    disp(k);
    devs_model2B=nan(length(c_0_test),length(tryD));
    devs_model2B_F2=nan(length(c_0_test),length(tryD));
    for i=1:length(c_0_test)
        for j=1:length(tryD)
            c_0=c_0_test(i);
            devs_model2B(i,j)=(V_F1(k)-c_0.*F_CDF(k).*(1-tryD(j)).*d_F1(k));
            devs_model2B_F2(i,j)=(V_F2(k)-c_0.*F_CDF(k).*d_F2(k)-c_0.*F_CDF(k).*tryD(j).*d_F1(k));
        end
    end
    atEachFreq{k}=min(devs_model2B,[],2)+min(devs_model2B_F2,[],2);
%     atEachFreq{k}=min(devs_model2B,[],2);
end
allErrs=zeros(size(atEachFreq{1}));
for k=1:length(freqs)-3
    allErrs=allErrs+atEachFreq{k};
end

% allErrs(~canUse)=nan;
[mi,mind]=nanmin(allErrs);
c_0=c_0_test(mind);
c_0=mean(V_F1)./(0.2.*mean(d_F1.*F_CDF));
% c_0=mean([c_0_test_start c_0_test_end]);
% Model 2B - filtering by CDF
D_fromeq1_2B=1-(V_F1./(c_0.*F_CDF.*d_F1)); % D as a function of f
D_fromeq2_2B=(V_F2-c_0.*F_CDF.*d_F2)./(c_0.*F_CDF.*d_F1); % D as a function of f
% V_F1_pred_2B=c_0.*F_CDF.*(1-D_fromeq2_2B).*d_F1;

% figure();
% plot(freqs,V_F1_pred,'Color','g');
% hold on;
% plot(freqs,V_F1,'Color','k');
% title('Model 2A no LP filtering');
% 
% figure();
% plot(freqs,V_F1_pred_2B,'Color','g');
% hold on;
% plot(freqs,V_F1,'Color','k');
% title('Model 2B w LP filtering');

figure();
semilogx(freqs,-D_fromeq1,'-k');
hold on;
semilogx(freqs,D_fromeq2,'-r');
title('Model 2A no LP filtering');

figure();
semilogx(freqs,D_fromeq1_2B,'-k');
hold on;
semilogx(freqs,D_fromeq2_2B,'-r');
title('Model 2B w LP filtering');

