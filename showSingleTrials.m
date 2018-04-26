

scrsz=get(0,'ScreenSize');
%figure('Position',[100 scrsz(4)/3.2-300 scrsz(3)*0.57 scrsz(4)/8]); 
%figure('Position',[0 scrsz(4)/3.2-600 scrsz(3)*0.57 scrsz(4)]); 
figure('Position',[0 scrsz(4)/3.2-600 scrsz(3)*0.57 scrsz(4)*(8/16)]); 

n=151;

%blue=[4 5 7 9 11 13 15 17 19 21 23 25 29 31 33 35 37 39 41 43 45 47]; daq33 
%blue=[3 5 7 10 11 13 15 17 19 21 23 25 27 29 33 35 37 39 41 43 45 47];

for i=(n-1)*5+1:(n-1)*5+1+4
    subplot(5,1,i-(n-1)*5);
    i-(n-1)*5
    if ledForSweeps(i)>0
    %if ismember(i,blue)
        c='b';
    else
        c='k';
    end
    plot(0:1/LFP_Fs:(1/LFP_Fs)*(size(LFPbySweep,2)-1),LFPbySweep(i,:),'Color',c,'LineWidth',1); 
    axis tight;
end