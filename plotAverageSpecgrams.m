function plotAverageSpecgrams(datadir)

placeStimOnsetAt=4; % seconds from trial onset
cutToMaxTime=12; % drop times above this in seconds
cutToMinTime=2; % drop times below this in seconds

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
    
    all_noTheta_noLED_Ntsr1.t=[];
    all_noTheta_noLED_Ntsr1.f=[];
    all_noTheta_noLED_Ntsr1.S=[];
    all_noTheta_LED_Ntsr1.t=[];
    all_noTheta_LED_Ntsr1.f=[];
    all_noTheta_LED_Ntsr1.S=[];
    all_theta_noLED_Ntsr1.t=[];
    all_theta_noLED_Ntsr1.f=[];
    all_theta_noLED_Ntsr1.S=[];
    all_theta_LED_Ntsr1.t=[];
    all_theta_LED_Ntsr1.f=[];
    all_theta_LED_Ntsr1.S=[];
    
    all_noTheta_noLED.HFa=[];
    all_noTheta_noLED.LFa=[];
    all_noTheta_noLED.F1amp=[];
    all_noTheta_LED.HFa=[];
    all_noTheta_LED.LFa=[];
    all_noTheta_LED.F1amp=[];
    all_theta_noLED.HFa=[];
    all_theta_noLED.LFa=[];
    all_theta_noLED.F1amp=[];
    all_theta_LED.HFa=[];
    all_theta_LED.LFa=[];
    all_theta_LED.F1amp=[];
    
    all_noTheta_noLED_Ntsr1.HFa=[];
    all_noTheta_noLED_Ntsr1.LFa=[];
    all_noTheta_noLED_Ntsr1.F1amp=[];
    all_noTheta_LED_Ntsr1.HFa=[];
    all_noTheta_LED_Ntsr1.LFa=[];
    all_noTheta_LED_Ntsr1.F1amp=[];
    all_theta_noLED_Ntsr1.HFa=[];
    all_theta_noLED_Ntsr1.LFa=[];
    all_theta_noLED_Ntsr1.F1amp=[];
    all_theta_LED_Ntsr1.HFa=[];
    all_theta_LED_Ntsr1.LFa=[];
    all_theta_LED_Ntsr1.F1amp=[];
    
    runningTallyNtsr1=0;
    runningTallyOther=0;
    for i=1:length(datadir)
        d=datadir{i};
        a=load([d '\' 'output.mat']);
        ntsr=a.output;
        a=load([d '\' 'wvfms.mat']);
        wvfms=a.wvfms;
        a=load([d '\' 'deepenough.mat']);
        deepenough=a.deepenough;
        if length(deepenough)==1
            deepenough=ones(1,length(wvfms.hw)).*deepenough;
        end
        a=load([d '\' 'stimWindow.mat']);
        stimWindow=a.stimWindow;
        
        a=load([d '\' 'noTheta_noLED']);
        noTheta_noLED=a.noTheta;
        a=load([d '\' 'noTheta_LED']);
        noTheta_LED=a.noTheta;
        a=load([d '\' 'theta_noLED']);
        theta_noLED=a.theta;
        a=load([d '\' 'theta_LED']);
        theta_LED=a.theta;
        
        [noTheta_noLED,noTheta_LED,theta_noLED,theta_LED]=alignToStimOnset(noTheta_noLED,noTheta_LED,theta_noLED,theta_LED,stimWindow,placeStimOnsetAt);
        
        if exist([d '\' 'classifyAsNtsr1.mat'],'file')
            classifyAlreadyExists=true;
            a=load([d '\' 'classifyAsNtsr1.mat']);
            classifyAsNtsr1=a.classifyAsNtsr1;
        else
            classifyAlreadyExists=false;
            wvfms.userIsFS=nan(1,length(wvfms.hw));
            classifyAsNtsr1=nan(1,length(wvfms.hw));
            for j=1:length(wvfms.hw)
                currhw=wvfms.hw(j);
                close all;
                figure();
                if length(size(wvfms.spikewvfms))==3
                    waveform=wvfms.spikewvfms(j,:,:);
                    waveform=waveform(1:end);
                    set(gcf,'Position',[100 300 800 500]);
                else
                    waveform=wvfms.spikewvfms(j,:);
                    set(gcf,'Position',[100 300 200 500]);
                end
                plot(waveform);
                if currhw<0.32*10^-3
                    % may be FS
                    answer=questdlg('Fast-spiking? Think so','Find RS vs FS','yes','no','yes');
                else
                    % may be RS
                    answer=questdlg('Fast-spiking? Doubt it','Find RS vs FS','yes','no','no');
                end
                if strcmp(answer,'yes')
                    wvfms.userIsFS(j)=1;
                elseif strcmp(answer,'no')
                    wvfms.userIsFS(j)=0;
                else
                    error('unrecognized button press');
                end
                if ntsr.putative_ntsr1_activity_pattern(j)==1 && wvfms.userIsFS(j)==0 && deepenough(j)==1
                    classifyAsNtsr1(j)=1;
                else
                    classifyAsNtsr1(j)=0;
                end
            end
        end
        if classifyAlreadyExists==false
            save([d '\' 'classifyAsNtsr1.mat'],'classifyAsNtsr1');
            isFs=wvfms.userIsFS;
            save([d '\' 'isFs.mat'],'isFs');
        end
        for j=1:length(wvfms.hw)
            if (j==1 && i==1) || (classifyAsNtsr1(j)==1 && isempty(all_noTheta_noLED_Ntsr1.S)) || (classifyAsNtsr1(j)==0 && isempty(all_noTheta_noLED.S))
                if classifyAsNtsr1(j)==1
                    all_noTheta_noLED_Ntsr1.t=noTheta_noLED.allS.t;
                    all_noTheta_noLED_Ntsr1.f=noTheta_noLED.allS.f;
                    all_noTheta_noLED_Ntsr1.S=noTheta_noLED.allS.S{j};
                    all_noTheta_LED_Ntsr1.t=noTheta_LED.allS.t;
                    all_noTheta_LED_Ntsr1.f=noTheta_LED.allS.f;
                    all_noTheta_LED_Ntsr1.S=noTheta_LED.allS.S{j};
                    all_theta_noLED_Ntsr1.t=theta_noLED.allS.t;
                    all_theta_noLED_Ntsr1.f=theta_noLED.allS.f;
                    all_theta_noLED_Ntsr1.S=theta_noLED.allS.S{j};
                    all_theta_LED_Ntsr1.t=theta_LED.allS.t;
                    all_theta_LED_Ntsr1.f=theta_LED.allS.f;
                    all_theta_LED_Ntsr1.S=theta_LED.allS.S{j};
                    runningTallyNtsr1=runningTallyNtsr1+1;
                    
                    all_noTheta_noLED_Ntsr1.HFa=noTheta_noLED.HFa;
                    all_noTheta_noLED_Ntsr1.LFa=noTheta_noLED.LFa;
                    all_noTheta_noLED_Ntsr1.F1amp=noTheta_noLED.F1amp;
                    all_noTheta_LED_Ntsr1.HFa=noTheta_LED.HFa;
                    all_noTheta_LED_Ntsr1.LFa=noTheta_LED.LFa;
                    all_noTheta_LED_Ntsr1.F1amp=noTheta_LED.F1amp;
                    all_theta_noLED_Ntsr1.HFa=theta_noLED.HFa;
                    all_theta_noLED_Ntsr1.LFa=theta_noLED.LFa;
                    all_theta_noLED_Ntsr1.F1amp=theta_noLED.F1amp;
                    all_theta_LED_Ntsr1.HFa=theta_LED.HFa;
                    all_theta_LED_Ntsr1.LFa=theta_LED.LFa;
                    all_theta_LED_Ntsr1.F1amp=theta_LED.F1amp;
                else
                    all_noTheta_noLED.t=noTheta_noLED.allS.t;
                    all_noTheta_noLED.f=noTheta_noLED.allS.f;
                    all_noTheta_noLED.S=noTheta_noLED.allS.S{j};
                    all_noTheta_LED.t=noTheta_LED.allS.t;
                    all_noTheta_LED.f=noTheta_LED.allS.f;
                    all_noTheta_LED.S=noTheta_LED.allS.S{j};
                    all_theta_noLED.t=theta_noLED.allS.t;
                    all_theta_noLED.f=theta_noLED.allS.f;
                    all_theta_noLED.S=theta_noLED.allS.S{j};
                    all_theta_LED.t=theta_LED.allS.t;
                    all_theta_LED.f=theta_LED.allS.f;
                    all_theta_LED.S=theta_LED.allS.S{j};
                    runningTallyOther=runningTallyOther+1;
                    
                    all_noTheta_noLED.HFa=noTheta_noLED.HFa;
                    all_noTheta_noLED.LFa=noTheta_noLED.LFa;
                    all_noTheta_noLED.F1amp=noTheta_noLED.F1amp;
                    all_noTheta_LED.HFa=noTheta_LED.HFa;
                    all_noTheta_LED.LFa=noTheta_LED.LFa;
                    all_noTheta_LED.F1amp=noTheta_LED.F1amp;
                    all_theta_noLED.HFa=theta_noLED.HFa;
                    all_theta_noLED.LFa=theta_noLED.LFa;
                    all_theta_noLED.F1amp=theta_noLED.F1amp;
                    all_theta_LED.HFa=theta_LED.HFa;
                    all_theta_LED.LFa=theta_LED.LFa;
                    all_theta_LED.F1amp=theta_LED.F1amp;
                end
            else
                if classifyAsNtsr1(j)==1
                    [noTheta_noLED.allS.t,noTheta_noLED.allS.f,noTheta_noLED.allS.S{j}]=cutToSize(all_noTheta_noLED_Ntsr1,noTheta_noLED.allS.t,noTheta_noLED.allS.f,noTheta_noLED.allS.S{j});
                    [noTheta_LED.allS.t,noTheta_LED.allS.f,noTheta_LED.allS.S{j}]=cutToSize(all_noTheta_LED_Ntsr1,noTheta_LED.allS.t,noTheta_LED.allS.f,noTheta_LED.allS.S{j});
                    [theta_noLED.allS.t,theta_noLED.allS.f,theta_noLED.allS.S{j}]=cutToSize(all_theta_noLED_Ntsr1,theta_noLED.allS.t,theta_noLED.allS.f,theta_noLED.allS.S{j});
                    [theta_LED.allS.t,theta_LED.allS.f,theta_LED.allS.S{j}]=cutToSize(all_theta_LED_Ntsr1,theta_LED.allS.t,theta_LED.allS.f,theta_LED.allS.S{j});
                    
                    tmp=cat(3,all_noTheta_noLED_Ntsr1.S,noTheta_noLED.allS.S{j}); 
                    all_noTheta_noLED_Ntsr1.S=nansum(tmp,3);
                    tmp=cat(3,all_noTheta_LED_Ntsr1.S,noTheta_LED.allS.S{j}); 
                    all_noTheta_LED_Ntsr1.S=nansum(tmp,3);
                    tmp=cat(3,all_theta_noLED_Ntsr1.S,theta_noLED.allS.S{j}); 
                    all_theta_noLED_Ntsr1.S=nansum(tmp,3);
                    tmp=cat(3,all_theta_LED_Ntsr1.S,theta_LED.allS.S{j}); 
                    all_theta_LED_Ntsr1.S=nansum(tmp,3);
                    runningTallyNtsr1=runningTallyNtsr1+1;
                    
                    if size(noTheta_noLED.HFa,2)>size(all_noTheta_noLED_Ntsr1.HFa,2)
                        all_noTheta_noLED_Ntsr1.HFa=[all_noTheta_noLED_Ntsr1.HFa; noTheta_noLED.HFa(:,1:size(all_noTheta_noLED_Ntsr1.HFa,2))];
                        all_noTheta_LED_Ntsr1.HFa=[all_noTheta_LED_Ntsr1.HFa; noTheta_LED.HFa(:,1:size(all_noTheta_LED_Ntsr1.HFa,2))];
                        all_theta_noLED_Ntsr1.HFa=[all_theta_noLED_Ntsr1.HFa; theta_noLED.HFa(:,1:size(all_theta_noLED_Ntsr1.HFa,2))];
                        all_theta_LED_Ntsr1.HFa=[all_theta_LED_Ntsr1.HFa; theta_LED.HFa(:,1:size(all_theta_LED_Ntsr1.HFa,2))];
                        
                        all_noTheta_noLED_Ntsr1.LFa=[all_noTheta_noLED_Ntsr1.LFa; noTheta_noLED.LFa(:,1:size(all_noTheta_noLED_Ntsr1.LFa,2))];
                        all_noTheta_LED_Ntsr1.LFa=[all_noTheta_LED_Ntsr1.LFa; noTheta_LED.LFa(:,1:size(all_noTheta_LED_Ntsr1.LFa,2))];
                        all_theta_noLED_Ntsr1.LFa=[all_theta_noLED_Ntsr1.LFa; theta_noLED.LFa(:,1:size(all_theta_noLED_Ntsr1.LFa,2))];
                        all_theta_LED_Ntsr1.LFa=[all_theta_LED_Ntsr1.LFa; theta_LED.LFa(:,1:size(all_theta_LED_Ntsr1.LFa,2))];
                    elseif size(noTheta_noLED.HFa,2)<size(all_noTheta_noLED_Ntsr1.HFa,2)
                        all_noTheta_noLED_Ntsr1.HFa=[all_noTheta_noLED_Ntsr1.HFa; [noTheta_noLED.HFa nan(size(noTheta_noLED.HFa,1),size(all_noTheta_noLED_Ntsr1.HFa,2)-size(noTheta_noLED.HFa,2))]];
                        all_noTheta_LED_Ntsr1.HFa=[all_noTheta_LED_Ntsr1.HFa; [noTheta_LED.HFa nan(size(noTheta_LED.HFa,1),size(all_noTheta_LED_Ntsr1.HFa,2)-size(noTheta_LED.HFa,2))]];
                        all_theta_noLED_Ntsr1.HFa=[all_theta_noLED_Ntsr1.HFa; [theta_noLED.HFa nan(size(theta_noLED.HFa,1),size(all_theta_noLED_Ntsr1.HFa,2)-size(theta_noLED.HFa,2))]];
                        all_theta_LED_Ntsr1.HFa=[all_theta_LED_Ntsr1.HFa; [theta_LED.HFa nan(size(theta_LED.HFa,1),size(all_theta_LED_Ntsr1.HFa,2)-size(theta_LED.HFa,2))]];
                        
                        all_noTheta_noLED_Ntsr1.LFa=[all_noTheta_noLED_Ntsr1.LFa; [noTheta_noLED.LFa nan(size(noTheta_noLED.LFa,1),size(all_noTheta_noLED_Ntsr1.LFa,2)-size(noTheta_noLED.LFa,2))]];
                        all_noTheta_LED_Ntsr1.LFa=[all_noTheta_LED_Ntsr1.LFa; [noTheta_LED.LFa nan(size(noTheta_LED.LFa,1),size(all_noTheta_LED_Ntsr1.LFa,2)-size(noTheta_LED.LFa,2))]];
                        all_theta_noLED_Ntsr1.LFa=[all_theta_noLED_Ntsr1.LFa; [theta_noLED.LFa nan(size(theta_noLED.LFa,1),size(all_theta_noLED_Ntsr1.LFa,2)-size(theta_noLED.LFa,2))]];
                        all_theta_LED_Ntsr1.LFa=[all_theta_LED_Ntsr1.LFa; [theta_LED.LFa nan(size(theta_LED.LFa,1),size(all_theta_LED_Ntsr1.LFa,2)-size(theta_LED.LFa,2))]];
                    end   
                else
                    [noTheta_noLED.allS.t,noTheta_noLED.allS.f,noTheta_noLED.allS.S{j}]=cutToSize(all_noTheta_noLED,noTheta_noLED.allS.t,noTheta_noLED.allS.f,noTheta_noLED.allS.S{j});
                    [noTheta_LED.allS.t,noTheta_LED.allS.f,noTheta_LED.allS.S{j}]=cutToSize(all_noTheta_LED,noTheta_LED.allS.t,noTheta_LED.allS.f,noTheta_LED.allS.S{j});
                    [theta_noLED.allS.t,theta_noLED.allS.f,theta_noLED.allS.S{j}]=cutToSize(all_theta_noLED,theta_noLED.allS.t,theta_noLED.allS.f,theta_noLED.allS.S{j});
                    [theta_LED.allS.t,theta_LED.allS.f,theta_LED.allS.S{j}]=cutToSize(all_theta_LED,theta_LED.allS.t,theta_LED.allS.f,theta_LED.allS.S{j});
                    
                    tmp=cat(3,all_noTheta_noLED.S,noTheta_noLED.allS.S{j}); 
                    all_noTheta_noLED.S=nansum(tmp,3);
                    tmp=cat(3,all_noTheta_LED.S,noTheta_LED.allS.S{j}); 
                    all_noTheta_LED.S=nansum(tmp,3);
                    tmp=cat(3,all_theta_noLED.S,theta_noLED.allS.S{j}); 
                    all_theta_noLED.S=nansum(tmp,3);
                    tmp=cat(3,all_theta_LED.S,theta_LED.allS.S{j}); 
                    all_theta_LED.S=nansum(tmp,3);
                    runningTallyOther=runningTallyOther+1;
                    
                    if size(noTheta_noLED.HFa,2)>size(all_noTheta_noLED.HFa,2)
                        all_noTheta_noLED.HFa=[all_noTheta_noLED.HFa; noTheta_noLED.HFa(:,1:size(all_noTheta_noLED.HFa,2))];
                        all_noTheta_LED.HFa=[all_noTheta_LED.HFa; noTheta_LED.HFa(:,1:size(all_noTheta_LED.HFa,2))];
                        all_theta_noLED.HFa=[all_theta_noLED.HFa; theta_noLED.HFa(:,1:size(all_theta_noLED.HFa,2))];
                        all_theta_LED.HFa=[all_theta_LED.HFa; theta_LED.HFa(:,1:size(all_theta_LED.HFa,2))];
                        
                        all_noTheta_noLED.LFa=[all_noTheta_noLED.LFa; noTheta_noLED.LFa(:,1:size(all_noTheta_noLED.LFa,2))];
                        all_noTheta_LED.LFa=[all_noTheta_LED.LFa; noTheta_LED.LFa(:,1:size(all_noTheta_LED.LFa,2))];
                        all_theta_noLED.LFa=[all_theta_noLED.LFa; theta_noLED.LFa(:,1:size(all_theta_noLED.LFa,2))];
                        all_theta_LED.LFa=[all_theta_LED.LFa; theta_LED.LFa(:,1:size(all_theta_LED.LFa,2))];
                    elseif size(noTheta_noLED.HFa,2)<size(all_noTheta_noLED.HFa,2)
                        all_noTheta_noLED.HFa=[all_noTheta_noLED.HFa; [noTheta_noLED.HFa nan(size(noTheta_noLED.HFa,1),size(all_noTheta_noLED.HFa,2)-size(noTheta_noLED.HFa,2))]];
                        all_noTheta_LED.HFa=[all_noTheta_LED.HFa; [noTheta_LED.HFa nan(size(noTheta_LED.HFa,1),size(all_noTheta_LED.HFa,2)-size(noTheta_LED.HFa,2))]];
                        all_theta_noLED.HFa=[all_theta_noLED.HFa; [theta_noLED.HFa nan(size(theta_noLED.HFa,1),size(all_theta_noLED.HFa,2)-size(theta_noLED.HFa,2))]];
                        all_theta_LED.HFa=[all_theta_LED.HFa; [theta_LED.HFa nan(size(theta_LED.HFa,1),size(all_theta_LED.HFa,2)-size(theta_LED.HFa,2))]];
                        
                        all_noTheta_noLED.LFa=[all_noTheta_noLED.LFa; [noTheta_noLED.LFa nan(size(noTheta_noLED.LFa,1),size(all_noTheta_noLED.LFa,2)-size(noTheta_noLED.LFa,2))]];
                        all_noTheta_LED.LFa=[all_noTheta_LED.LFa; [noTheta_LED.LFa nan(size(noTheta_LED.LFa,1),size(all_noTheta_LED.LFa,2)-size(noTheta_LED.LFa,2))]];
                        all_theta_noLED.LFa=[all_theta_noLED.LFa; [theta_noLED.LFa nan(size(theta_noLED.LFa,1),size(all_theta_noLED.LFa,2)-size(theta_noLED.LFa,2))]];
                        all_theta_LED.LFa=[all_theta_LED.LFa; [theta_LED.LFa nan(size(theta_LED.LFa,1),size(all_theta_LED.LFa,2)-size(theta_LED.LFa,2))]];
                    end   
                end
            end
        end
    end
    all_noTheta_noLED_Ntsr1.S=all_noTheta_noLED_Ntsr1.S/runningTallyNtsr1;
    all_noTheta_LED_Ntsr1.S=all_noTheta_LED_Ntsr1.S/runningTallyNtsr1;
    all_theta_noLED_Ntsr1.S=all_theta_noLED_Ntsr1.S/runningTallyNtsr1;
    all_theta_LED_Ntsr1.S=all_theta_LED_Ntsr1.S/runningTallyNtsr1;
    all_noTheta_noLED.S=all_noTheta_noLED.S/runningTallyOther;
    all_noTheta_LED.S=all_noTheta_LED.S/runningTallyOther;
    all_theta_noLED.S=all_theta_noLED.S/runningTallyOther;
    all_theta_LED.S=all_theta_LED.S/runningTallyOther;
else
    disp('datadir should be cell array');
    return
end

figure(); 
imagesc(all_noTheta_noLED_Ntsr1.t(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_noTheta_noLED_Ntsr1.f(all_noTheta_noLED_Ntsr1.f<=50),all_noTheta_noLED_Ntsr1.S(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_noTheta_noLED_Ntsr1.f<=50)');
title('noTheta noLED Ntsr1');

figure(); 
temp=all_noTheta_noLED_Ntsr1.S(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_noTheta_noLED_Ntsr1.f<=50);
temp=temp./repmat(nansum(temp,2),1,size(temp,2));
imagesc(all_noTheta_noLED_Ntsr1.t(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_noTheta_noLED_Ntsr1.f(all_noTheta_noLED_Ntsr1.f<=50),temp');
title('noTheta noLED Ntsr1 normalized');

figure(); 
imagesc(all_noTheta_LED_Ntsr1.t(~isnan(all_noTheta_LED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_noTheta_LED_Ntsr1.f(all_noTheta_LED_Ntsr1.f<=50),all_noTheta_LED_Ntsr1.S(~isnan(all_noTheta_LED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_noTheta_LED_Ntsr1.f<=50)');
title('noTheta LED Ntsr1');

figure(); 
imagesc(all_theta_noLED_Ntsr1.t(~isnan(all_theta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_theta_noLED_Ntsr1.f(all_theta_noLED_Ntsr1.f<=50),all_theta_noLED_Ntsr1.S(~isnan(all_theta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_theta_noLED_Ntsr1.f<=50)');
title('theta noLED Ntsr1');

figure(); 
imagesc(all_theta_LED_Ntsr1.t(~isnan(all_theta_LED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_theta_LED_Ntsr1.f(all_theta_LED_Ntsr1.f<=50),all_theta_LED_Ntsr1.S(~isnan(all_theta_LED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_theta_LED_Ntsr1.f<=50)');
title('theta LED Ntsr1');


figure(); 
imagesc(all_noTheta_noLED.t(~isnan(all_noTheta_noLED.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_noTheta_noLED.f(all_noTheta_noLED.f<=50),all_noTheta_noLED.S(~isnan(all_noTheta_noLED.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_noTheta_noLED.f<=50)');
title('noTheta noLED');

figure(); 
imagesc(all_noTheta_LED.t(~isnan(all_noTheta_LED.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_noTheta_LED.f(all_noTheta_LED.f<=50),all_noTheta_LED.S(~isnan(all_noTheta_LED.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_noTheta_LED.f<=50)');
title('noTheta LED');

figure(); 
imagesc(all_theta_noLED.t(~isnan(all_theta_noLED.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_theta_noLED.f(all_theta_noLED.f<=50),all_theta_noLED.S(~isnan(all_theta_noLED.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_theta_noLED.f<=50)');
title('theta noLED');

figure(); 
imagesc(all_theta_LED.t(~isnan(all_theta_LED.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_theta_LED.f(all_theta_LED.f<=50),all_theta_LED.S(~isnan(all_theta_LED.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_theta_LED.f<=50)');
title('theta LED');

plotWStderr(all_noTheta_noLED_Ntsr1.t,all_noTheta_noLED_Ntsr1.LFa,all_noTheta_noLED_Ntsr1.HFa,'k','r',1);
title('Ntsr1 No Theta No LED LF vs HF');

plotWStderr(all_theta_noLED_Ntsr1.t,all_theta_noLED_Ntsr1.LFa,all_theta_noLED_Ntsr1.HFa,'k','r',1);
title('Ntsr1 Theta No LED LF vs HF');

plotWStderr(all_noTheta_noLED.t,all_noTheta_noLED.LFa,all_noTheta_noLED.HFa,'k','r',1);
title('No Theta No LED LF vs HF');

plotWStderr(all_theta_noLED.t,all_theta_noLED.LFa,all_theta_noLED.HFa,'k','r',1);
title('Theta No LED LF vs HF');

plotWStderr(all_theta_noLED.t,all_noTheta_noLED_Ntsr1.LFa./all_noTheta_noLED_Ntsr1.HFa,all_noTheta_noLED.LFa./all_noTheta_noLED.HFa,'k','r',1);
title('No theta no LED, Ratio LF to HF -- Ntsr1 vs Non-Ntsr1');

plotWStderr(all_theta_noLED.t,all_theta_noLED_Ntsr1.LFa./all_theta_noLED_Ntsr1.HFa,all_theta_noLED.LFa./all_theta_noLED.HFa,'k','r',1);
title('Theta no LED, Ratio LF to HF -- Ntsr1 vs Non-Ntsr1');

figure(); 
temp=all_noTheta_noLED_Ntsr1.S(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_noTheta_noLED_Ntsr1.f<=50);
temp=temp./repmat(nansum(temp,2),1,size(temp,2));
imagesc(all_noTheta_noLED_Ntsr1.t(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_noTheta_noLED_Ntsr1.f(all_noTheta_noLED_Ntsr1.f<=50),temp');
title('noTheta noLED Ntsr1 normalized');

end

function plotWStderr(x,y1,y2,c1,c2,doFill)

figure();
plot(x,nanmean(y1,1),'Color',c1);
hold on;
if doFill==1
    fill([x fliplr(x)],[nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)) fliplr(nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)))],[0.5 0.5 0.5]);
end
plot(x,nanmean(y1,1),'Color',c1);
plot(x,nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);
plot(x,nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);

plot(x,nanmean(y2,1),'Color',c2);
if doFill==1
    fill([x fliplr(x)],[nanmean(y2,1)+nanstd(y2,[],1)./sqrt(size(y2,1)) fliplr(nanmean(y2,1)-nanstd(y2,[],1)./sqrt(size(y2,1)))],[0.1 0.7 0.5]);
end
plot(x,nanmean(y2,1),'Color',c2);
plot(x,nanmean(y2,1)+nanstd(y2,[],1)./sqrt(size(y2,1)),'Color',c2);
plot(x,nanmean(y2,1)-nanstd(y2,[],1)./sqrt(size(y2,1)),'Color',c2);

plot(x,nanmean(y1,1),'Color',c1);
legend({'LFa','','','','','HFa','','','',''});
    
end

function [newData_t,newData_f,newData_S]=cutToSize(fitTo,newData_t,newData_f,newData_S)

t_size=length(fitTo.t);
f_size=length(fitTo.f);
S_size=size(fitTo.S);

newData_S=newData_S(newData_t>=nanmin(fitTo.t) & newData_t<=nanmax(fitTo.t),newData_f>=nanmin(fitTo.f) & newData_f<=nanmax(fitTo.f));
newData_t=newData_t(newData_t>=nanmin(fitTo.t) & newData_t<=nanmax(fitTo.t));
newData_f=newData_f(newData_f>=nanmin(fitTo.f) & newData_f<=nanmax(fitTo.f));

if length(newData_t)>t_size
    newData_S=newData_S(1:t_size,:);
    newData_t=newData_t(1:t_size);
elseif length(newData_t)<t_size
    newData_S=[newData_S; nan(t_size-length(newData_t),size(newData_S,2))];
    newData_t=[newData_t nan(1,t_size-length(newData_t))];
end
if length(newData_f)>f_size
    newData_S=newData_S(:,1:f_size);
    newData_f=newData_f(1:f_size);
elseif length(newData_f)<f_size
    newData_S=[newData_S; nan(size(newData_S,1),f_size-length(newData_f))];
    newData_f=[newData_f nan(1,f_size-length(newData_f))];
end

end

function [noTheta_noLED,noTheta_LED,theta_noLED,theta_LED]=alignToStimOnset(noTheta_noLED,noTheta_LED,theta_noLED,theta_LED,stimWindow,placeStimOnsetAt)

stimOnset=stimWindow(1);
if abs(stimOnset-placeStimOnsetAt)<0.01
    return
end

noTheta_noLED=alignEach(noTheta_noLED,stimOnset,placeStimOnsetAt);
noTheta_LED=alignEach(noTheta_LED,stimOnset,placeStimOnsetAt);
theta_noLED=alignEach(theta_noLED,stimOnset,placeStimOnsetAt);
theta_LED=alignEach(theta_LED,stimOnset,placeStimOnsetAt);

end

function dataStruct=alignEach(dataStruct,stimOnset,placeStimOnsetAt)

% assume same binning for all
shiftByInds=floor((stimOnset-placeStimOnsetAt)/(dataStruct.allS.t(2)-dataStruct.allS.t(1)));

if shiftByInds>0
    % shift forward by shiftByInds
    for i=1:length(dataStruct.allS.S)
        temp=dataStruct.allS.S{i};
        temp=[temp(shiftByInds:end,:); nan(shiftByInds-1,size(temp,2))];
        dataStruct.allS.S{i}=temp;
    end
    dataStruct.HFa=[dataStruct.HFa(:,shiftByInds:end) nan(size(dataStruct.HFa,1),shiftByInds-1)];
    dataStruct.LFa=[dataStruct.LFa(:,shiftByInds:end) nan(size(dataStruct.LFa,1),shiftByInds-1)];
    dataStruct.F1amp=[dataStruct.F1amp(:,shiftByInds:end) nan(size(dataStruct.F1amp,1),shiftByInds-1)];
    dataStruct.allpower=[dataStruct.allpower(:,shiftByInds:end) nan(size(dataStruct.allpower,1),shiftByInds-1)];
elseif shiftByInds<0
    % shift backward by shiftByInds
    for i=1:length(dataStruct.allS.S)
        temp=dataStruct.allS.S{i};
        temp=[nan(shiftByInds,size(temp,2)); temp(1:end-shiftByInds,:)];
        dataStruct.allS.S{i}=temp;
    end
    dataStruct.HFa=[nan(size(dataStruct.HFa,1),shiftByInds) dataStruct.HFa(:,1:end-shiftByInds)];
    dataStruct.LFa=[nan(size(dataStruct.LFa,1),shiftByInds) dataStruct.LFa(:,1:end-shiftByInds)];
    dataStruct.F1amp=[nan(size(dataStruct.F1amp,1),shiftByInds) dataStruct.F1amp(:,1:end-shiftByInds)];
    dataStruct.allpower=[nan(size(dataStruct.allpower,1),shiftByInds) dataStruct.allpower(:,1:end-shiftByInds)];
end

end