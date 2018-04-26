function getF1effects_bySpatFreq(filedir_ind)

% Spat Freq for each expt
dd{1}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake310\Lower Theta Thresh';
dd{2}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake311';
dd{3}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake312\Lower Theta Thresh';
dd{4}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake313';
dd{5}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake314\Lower Theta Thresh';
dd{6}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake315';
dd{7}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake316\Lower Theta Thresh';
dd{8}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake317\Lower Theta Thresh';
dd{9}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake318';
dd{10}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake319';
dd{11}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake320';
dd{12}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake321';
dd{13}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake322';
dd{14}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake327\Lower Theta Thresh';
dd{15}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake328\Lower Theta Thresh';
dd{16}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake329\Lower Theta Thresh';
dd{17}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake330\Lower Theta Thresh';
dd{18}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake373\Lower Theta Thresh';
dd{19}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake376\Lower Theta Thresh';
dd{20}='C:\Users\Admin\Documents\Figures\Figure 1c\Mawake377\Lower Theta Thresh';

sf{1}=[0.03,0.03,0.03];
sf{2}=[0.03,0.03,0.03,0.03];
sf{3}=[0.03,0.03,0.03];
sf{4}=[0.03,0.03,0.03];
sf{5}=[0.03,0.03,0.03,0.03];
sf{6}=[0.03,0.03,0.03];
sf{7}=[0.02,0.08,0.04];
sf{8}=[0.02,0.08,0.04];
sf{9}=[0.02,0.08,0.04];
sf{10}=[0.04,0.08,0.02];
sf{11}=[0.03,0.08];
sf{12}=[0.08,0.03];
sf{13}=[0.03,0.08];
sf{14}=[0.02,0.04];
sf{15}=[0.02,0.04];
sf{16}=[0.02,0.04];
sf{17}=[0.03,0.06];
sf{18}=[0.02,0.04,0.08];
sf{19}=[0.02,0.04,0.08];
sf{20}=[0.02,0.04,0.08];


filedir=dd{filedir_ind};
currsf=sf{filedir_ind};

a=load([filedir '\noTheta_trialAv_temp_noLED.mat']);
noTheta_noLED=a.noTheta_trialAv_temp;
a=load([filedir '\theta_trialAv_temp_noLED.mat']);
theta_noLED=a.theta_trialAv_temp;
a=load([filedir '\noTheta_trialAv_temp_LED.mat']);
noTheta_LED=a.noTheta_trialAv_temp;
a=load([filedir '\theta_trialAv_temp_LED.mat']);
theta_LED=a.theta_trialAv_temp;

if length(unique(currsf))==1 && length(currsf)>1
    % Average
    temp_noTheta_noLED=zeros(size(noTheta_noLED(1).F1amp));
    temp_theta_noLED=zeros(size(theta_noLED(1).F1amp));
    temp_noTheta_LED=zeros(size(noTheta_LED(1).F1amp));
    temp_theta_LED=zeros(size(theta_LED(1).F1amp));
    for i=1:length(currsf)
        tempie(:,:,1)=temp_noTheta_noLED;
        tempie(:,:,2)=noTheta_noLED(i).F1amp;
        temp_noTheta_noLED=nansum(tempie,3);
%         temp_noTheta_noLED=temp_noTheta_noLED+noTheta_noLED(i).F1amp;
        tempie(:,:,1)=temp_theta_noLED;
        tempie(:,:,2)=theta_noLED(i).F1amp;
        temp_theta_noLED=nansum(tempie,3);
%         temp_theta_noLED=temp_theta_noLED+theta_noLED(i).F1amp;
        tempie(:,:,1)=temp_noTheta_LED;
        tempie(:,:,2)=noTheta_LED(i).F1amp;
        temp_noTheta_LED=nansum(tempie,3);
%         temp_noTheta_LED=temp_noTheta_LED+noTheta_LED(i).F1amp;
        tempie(:,:,1)=temp_theta_LED;
        tempie(:,:,2)=theta_LED(i).F1amp;
        temp_theta_LED=nansum(tempie,3);
%         temp_theta_LED=temp_theta_LED+theta_LED(i).F1amp;
    end
    temp_noTheta_noLED=temp_noTheta_noLED./length(currsf);
    temp_theta_noLED=temp_theta_noLED./length(currsf);
    temp_noTheta_LED=temp_noTheta_LED./length(currsf);
    temp_theta_LED=temp_theta_LED./length(currsf);
    
    noTheta_noLED(1).F1amp=temp_noTheta_noLED;
    theta_noLED(1).F1amp=temp_theta_noLED;
    noTheta_LED(1).F1amp=temp_noTheta_LED;
    theta_LED(1).F1amp=temp_theta_LED;
    
    currdir=[filedir '\SF0' num2str(currsf(1)*100)];
    mkdir(currdir);
    noTheta_trialAv=noTheta_noLED(1);
    theta_trialAv=theta_noLED(1);
    save([currdir '\noTheta_trialAv_noLED.mat'],'noTheta_trialAv');
    save([currdir '\theta_trialAv_noLED.mat'],'theta_trialAv');
    noTheta_trialAv=noTheta_LED(1);
    theta_trialAv=theta_LED(1);
    save([currdir '\noTheta_trialAv_LED.mat'],'noTheta_trialAv');
    save([currdir '\theta_trialAv_LED.mat'],'theta_trialAv');
else
    for i=1:length(currsf)
        noTheta_trialAv_noLED=noTheta_noLED(i);
        theta_trialAv_noLED=theta_noLED(i);
        noTheta_trialAv_LED=noTheta_LED(i);
        theta_trialAv_LED=theta_LED(i);
        currdir=[filedir '\SF0' num2str(currsf(i)*100)];
        mkdir(currdir);
        noTheta_trialAv=noTheta_noLED(i);
        theta_trialAv=theta_noLED(i);
        save([currdir '\noTheta_trialAv_noLED.mat'],'noTheta_trialAv');
        save([currdir '\theta_trialAv_noLED.mat'],'theta_trialAv');
        noTheta_trialAv=noTheta_LED(i);
        theta_trialAv=theta_LED(i);
        save([currdir '\noTheta_trialAv_LED.mat'],'noTheta_trialAv');
        save([currdir '\theta_trialAv_LED.mat'],'theta_trialAv');
    end
end