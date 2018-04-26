function plotThetaNoThetaSpecgrams(datadir)

if iscell(datadir)
    all_noTheta_trialAv_noLED.t=[];
    all_noTheta_trialAv_noLED.f=[];
    all_noTheta_trialAv_noLED.S=[];
    all_noTheta_trialAv_LED.t=[];
    all_noTheta_trialAv_LED.f=[];
    all_noTheta_trialAv_LED.S=[];
    all_theta_trialAv_noLED.t=[];
    all_theta_trialAv_noLED.f=[];
    all_theta_trialAv_noLED.S=[];
    all_theta_trialAv_LED.t=[];
    all_theta_trialAv_LED.f=[];
    all_theta_trialAv_LED.S=[];
    for i=1:length(datadir)
        d=datadir{i};
        a=load([d '\' 'noTheta_trialAv_temp_noLED']);
        noTheta_trialAv_noLED=a.noTheta_trialAv_temp;
        avSpecgram=avS(noTheta_trialAv_noLED);
        if i==1
            all_noTheta_trialAv_noLED.t=avSpecgram.t;
            all_noTheta_trialAv_noLED.f=avSpecgram.f;
        end
        all_noTheta_trialAv_noLED.S=[all_noTheta_trialAv_noLED.S avSpecgram.S];
        
        a=load([d '\' 'noTheta_trialAv_temp_LED']);
        noTheta_trialAv_LED=a.noTheta_trialAv_temp;
        avSpecgram=avS(noTheta_trialAv_LED);
        if i==1
            all_noTheta_trialAv_LED.t=avSpecgram.t;
            all_noTheta_trialAv_LED.f=avSpecgram.f;
        end
        all_noTheta_trialAv_LED.S=[all_noTheta_trialAv_LED.S avSpecgram.S];
        
        a=load([d '\' 'theta_trialAv_temp_noLED']);
        theta_trialAv_noLED=a.theta_trialAv_temp;
        avSpecgram=avS(theta_trialAv_noLED);
        if i==1
            all_theta_trialAv_noLED.t=avSpecgram.t;
            all_theta_trialAv_noLED.f=avSpecgram.f;
        end
        all_theta_trialAv_noLED.S=[all_theta_trialAv_noLED.S avSpecgram.S];

        a=load([d '\' 'theta_trialAv_temp_LED']);
        theta_trialAv_LED=a.theta_trialAv_temp;
        avSpecgram=avS(theta_trialAv_LED);
        if i==1
            all_theta_trialAv_LED.t=avSpecgram.t;
            all_theta_trialAv_LED.f=avSpecgram.f;
        end
        all_theta_trialAv_LED.S=[all_theta_trialAv_LED.S avSpecgram.S];
    end
    noTheta_trialAv_noLED=all_noTheta_trialAv_noLED;
    noTheta_trialAv_LED=all_noTheta_trialAv_LED;
    theta_trialAv_noLED=all_theta_trialAv_noLED;
    theta_trialAv_LED=all_theta_trialAv_LED;
end

sumS=zeros(size(noTheta_trialAv_noLED.S{1}));
tally=0;
% allF1s_noTheta_noLED=nan(
for i=1:length(noTheta_trialAv_noLED.S)
    if ~isnan(noTheta_trialAv_noLED.S{i})
        sumS=sumS+noTheta_trialAv_noLED.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(noTheta_trialAv_noLED.t,noTheta_trialAv_noLED.f(noTheta_trialAv_noLED.f<=30),sumS(:,noTheta_trialAv_noLED.f<=30)');
title('noTheta_trialAv_noLED');

sumS=zeros(size(noTheta_trialAv_LED.S{1}));
tally=0;
for i=1:length(noTheta_trialAv_LED.S)
    if ~isnan(noTheta_trialAv_LED.S{i})
        sumS=sumS+noTheta_trialAv_LED.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(noTheta_trialAv_LED.t,noTheta_trialAv_LED.f(noTheta_trialAv_LED.f<=30),sumS(:,noTheta_trialAv_LED.f<=30)');
title('noTheta_trialAv_LED');

sumS=zeros(size(theta_trialAv_noLED.S{1}));
tally=0;
for i=1:length(theta_trialAv_noLED.S)
    if ~isnan(theta_trialAv_noLED.S{i})
        sumS=sumS+theta_trialAv_noLED.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(theta_trialAv_noLED.t,theta_trialAv_noLED.f(theta_trialAv_noLED.f<=30),sumS(:,theta_trialAv_noLED.f<=30)');
title('theta_trialAv_noLED');

sumS=zeros(size(theta_trialAv_LED.S{1}));
tally=0;
for i=1:length(theta_trialAv_LED.S)
    if ~isnan(theta_trialAv_LED.S{i})
        sumS=sumS+theta_trialAv_LED.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(theta_trialAv_LED.t,theta_trialAv_LED.f(theta_trialAv_LED.f<=30),sumS(:,theta_trialAv_LED.f<=30)');
title('theta_trialAv_LED');







if iscell(datadir)
    all_noTheta_noLED.t=[];
    all_noTheta_noLED.f=[];
    all_noTheta_noLED.S=[];
    all_noTheta_LED.t=[];
    all_noTheta_LED.f=[];
    all_noTheta_LED.S=[];
    all_theta_noLED.t=[];
    all_theta_noLED.f=[];
    all_theta_noLED.S=[];
    all_theta_LED.t=[];
    all_theta_LED.f=[];
    all_theta_LED.S=[];
    for i=1:length(datadir)
        d=datadir{i};
        a=load([d '\' 'noTheta_noLED']);
        noTheta_noLED=a.noTheta;
        avSpecgram=noTheta_noLED.allS;
        if i==1
            all_noTheta_noLED.t=avSpecgram.t;
            all_noTheta_noLED.f=avSpecgram.f;
        end
        all_noTheta_noLED.S=[all_noTheta_noLED.S avSpecgram.S];
        
        a=load([d '\' 'noTheta_LED']);
        noTheta_LED=a.noTheta;
        avSpecgram=noTheta_LED.allS;
        if i==1
            all_noTheta_LED.t=avSpecgram.t;
            all_noTheta_LED.f=avSpecgram.f;
        end
        all_noTheta_LED.S=[all_noTheta_LED.S avSpecgram.S];
        
        a=load([d '\' 'theta_noLED']);
        theta_noLED=a.theta;
        avSpecgram=theta_noLED.allS;
        if i==1
            all_theta_noLED.t=avSpecgram.t;
            all_theta_noLED.f=avSpecgram.f;
        end
        all_theta_noLED.S=[all_theta_noLED.S avSpecgram.S];

        a=load([d '\' 'theta_LED']);
        theta_LED=a.theta;
        avSpecgram=theta_LED.allS;
        if i==1
            all_theta_LED.t=avSpecgram.t;
            all_theta_LED.f=avSpecgram.f;
        end
        all_theta_LED.S=[all_theta_LED.S avSpecgram.S];
    end
    noTheta_noLED=all_noTheta_noLED;
    noTheta_LED=all_noTheta_LED;
    theta_noLED=all_theta_noLED;
    theta_LED=all_theta_LED;
end

sumS=zeros(size(noTheta_noLED.S{1}));
tally=0;
for i=1:length(noTheta_noLED.S)
    if ~isnan(noTheta_noLED.S{i})
        sumS=sumS+noTheta_noLED.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(noTheta_noLED.t,noTheta_noLED.f(noTheta_noLED.f<=30),sumS(:,noTheta_noLED.f<=30)');
title('noTheta_noLED');

sumS=zeros(size(noTheta_LED.S{1}));
tally=0;
for i=1:length(noTheta_LED.S)
    if ~isnan(noTheta_LED.S{i})
        sumS=sumS+noTheta_LED.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(noTheta_LED.t,noTheta_LED.f(noTheta_LED.f<=30),sumS(:,noTheta_LED.f<=30)');
title('noTheta_LED');

sumS=zeros(size(theta_noLED.S{1}));
tally=0;
for i=1:length(theta_noLED.S)
    if ~isnan(theta_noLED.S{i})
        sumS=sumS+theta_noLED.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(theta_noLED.t,theta_noLED.f(theta_noLED.f<=30),sumS(:,theta_noLED.f<=30)');
title('theta_noLED');

sumS=zeros(size(theta_LED.S{1}));
tally=0;
for i=1:length(theta_LED.S)
    if ~isnan(theta_LED.S{i})
        sumS=sumS+theta_LED.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(theta_LED.t,theta_LED.f(theta_LED.f<=30),sumS(:,theta_LED.f<=30)');
title('theta_LED');


end

function avSpecgram=avS(temp_struct)

avSpecgram.t=[];
avSpecgram.f=[];
avSpecgram.S=[];
if length(temp_struct)==1
    avSpecgram.t=temp_struct.allS.t;
    avSpecgram.f=temp_struct.allS.f;
    avSpecgram.S=temp_struct.allS.S;
else
    avSpecgram.t=temp_struct(1).allS.t;
    avSpecgram.f=temp_struct(1).allS.f;
    avSpecgram.S=cell(1,length(temp_struct(1).allS.S));
    for i=1:length(temp_struct(1).allS.S)
        avSpecgram.S{i}=zeros(size(temp_struct(1).allS.S{1}));
    end
    for i=1:length(temp_struct)
        for j=1:length(temp_struct(1).allS.S)
            avSpecgram.S{j}=avSpecgram.S{j}+temp_struct(i).allS.S{j};
        end
    end
    for i=1:length(avSpecgram.S)
        avSpecgram.S{i}=avSpecgram.S{i}./length(temp_struct);
    end
end

end


