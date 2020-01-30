function plotAverageSpecgrams(datadir)

placeStimOnsetAt=4; % seconds from trial onset
cutToMaxTime=12; % drop times above this in seconds
cutToMinTime=2; % drop times below this in seconds
normEachSpecgram=false; % normalize each spectrogram to its total power
takeFSForNonNtsr1=0;
spontWindow=[2 3.5];
evokedWindow=[4.5 6];

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
    all_noTheta_noLED.spont_spec=[];
    all_noTheta_noLED.evoked_spec=[];
    all_noTheta_LED.HFa=[];
    all_noTheta_LED.LFa=[];
    all_noTheta_LED.F1amp=[];
    all_noTheta_LED.spont_spec=[];
    all_noTheta_LED.evoked_spec=[];
    all_theta_noLED.HFa=[];
    all_theta_noLED.LFa=[];
    all_theta_noLED.F1amp=[];
    all_theta_noLED.spont_spec=[];
    all_theta_noLED.evoked_spec=[];
    all_theta_LED.HFa=[];
    all_theta_LED.LFa=[];
    all_theta_LED.F1amp=[];
    all_theta_LED.spont_spec=[];
    all_theta_LED.evoked_spec=[];
    
    all_noTheta_noLED_Ntsr1.HFa=[];
    all_noTheta_noLED_Ntsr1.LFa=[];
    all_noTheta_noLED_Ntsr1.F1amp=[];
    all_noTheta_noLED_Ntsr1.spont_spec=[];
    all_noTheta_noLED_Ntsr1.evoked_spec=[];
    all_noTheta_LED_Ntsr1.HFa=[];
    all_noTheta_LED_Ntsr1.LFa=[];
    all_noTheta_LED_Ntsr1.F1amp=[];
    all_noTheta_LED_Ntsr1.spont_spec=[];
    all_noTheta_LED_Ntsr1.evoked_spec=[];
    all_theta_noLED_Ntsr1.HFa=[];
    all_theta_noLED_Ntsr1.LFa=[];
    all_theta_noLED_Ntsr1.F1amp=[];
    all_theta_noLED_Ntsr1.spont_spec=[];
    all_theta_noLED_Ntsr1.evoked_spec=[];
    all_theta_LED_Ntsr1.HFa=[];
    all_theta_LED_Ntsr1.LFa=[];
    all_theta_LED_Ntsr1.F1amp=[];
    all_theta_LED_Ntsr1.spont_spec=[];
    all_theta_LED_Ntsr1.evoked_spec=[];
    
    all_noTheta_noLED_Ntsr1.psth=[];
    all_noTheta_LED_Ntsr1.psth=[];
    all_theta_noLED_Ntsr1.psth=[];
    all_theta_LED_Ntsr1.psth=[];
    all_noTheta_noLED.psth=[];
    all_noTheta_LED.psth=[];
    all_theta_noLED.psth=[];
    all_theta_LED.psth=[];
    psth_t=[];
    all_p_noTheta_noLED_Ntsr1=[];
    all_p_theta_noLED_Ntsr1=[];
    
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
        
        a=load([d '\' 'dLGNpsth']);
        psth=a.dLGNpsth;
        a=load([d '\' 'noThetaTrials']);
        noThetaTrials=a.noThetaTrials;
        
        [noTheta_noLED,noTheta_LED,theta_noLED,theta_LED,psth]=alignToStimOnset(noTheta_noLED,noTheta_LED,theta_noLED,theta_LED,stimWindow,placeStimOnsetAt,psth);
        
        if exist([d '\' 'classifyAsNtsr1.mat'],'file')
            classifyAlreadyExists=true;
            a=load([d '\' 'classifyAsNtsr1.mat']);
            classifyAsNtsr1=a.classifyAsNtsr1;
            a=load([d '\' 'isFs.mat']);
            isFs=a.isFs;
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
        if normEachSpecgram==true
            for j=1:length(wvfms.hw)
%                 noTheta_noLED.allS.S{j}=noTheta_noLED.allS.S{j}./nansum(nansum(noTheta_noLED.allS.S{j}));
%                 noTheta_LED.allS.S{j}=noTheta_LED.allS.S{j}./nansum(nansum(noTheta_LED.allS.S{j}));
%                 theta_noLED.allS.S{j}=theta_noLED.allS.S{j}./nansum(nansum(theta_noLED.allS.S{j}));
%                 theta_LED.allS.S{j}=theta_LED.allS.S{j}./nansum(nansum(theta_LED.allS.S{j}));
                
%                 if nansum(nansum(noTheta_noLED.allS.S{j}))<10*10*10^7
%                     noTheta_noLED.allS.S{j}=zeros(size(noTheta_noLED.allS.S{j}));
%                 end
%                 if nansum(nansum(noTheta_LED.allS.S{j}))<10*10*10^7
%                     noTheta_LED.allS.S{j}=zeros(size(noTheta_LED.allS.S{j}));
%                 end
%                 if nansum(nansum(theta_noLED.allS.S{j}))<10*10*10^7
%                     theta_noLED.allS.S{j}=zeros(size(theta_noLED.allS.S{j}));
%                 end
%                 if nansum(nansum(noTheta_noLED.allS.S{j}))<10*10*10^7
%                     theta_LED.allS.S{j}=zeros(size(theta_LED.allS.S{j}));
%                 end
                
                temp=noTheta_noLED.allS.S{j};
                temp=temp./repmat(nansum(temp,2),1,size(temp,2));
                noTheta_noLED.allS.S{j}=temp;
                temp=noTheta_LED.allS.S{j};
                temp=temp./repmat(nansum(temp,2),1,size(temp,2));
                noTheta_LED.allS.S{j}=temp;
                temp=theta_noLED.allS.S{j};
                temp=temp./repmat(nansum(temp,2),1,size(temp,2));
                theta_noLED.allS.S{j}=temp;
                temp=theta_LED.allS.S{j};
                temp=temp./repmat(nansum(temp,2),1,size(temp,2));
                theta_LED.allS.S{j}=temp;
            end
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
                    
                    all_noTheta_noLED_Ntsr1.spont_spec=grabSpecInTimeWindow(noTheta_noLED.allS.S{j},noTheta_noLED.allS.t,spontWindow);
                    all_noTheta_noLED_Ntsr1.evoked_spec=grabSpecInTimeWindow(noTheta_noLED.allS.S{j},noTheta_noLED.allS.t,evokedWindow);
                    all_noTheta_LED_Ntsr1.spont_spec=grabSpecInTimeWindow(noTheta_LED.allS.S{j},noTheta_LED.allS.t,spontWindow);
                    all_noTheta_LED_Ntsr1.evoked_spec=grabSpecInTimeWindow(noTheta_LED.allS.S{j},noTheta_LED.allS.t,evokedWindow);
                    all_theta_noLED_Ntsr1.spont_spec=grabSpecInTimeWindow(theta_noLED.allS.S{j},theta_noLED.allS.t,spontWindow);
                    all_theta_noLED_Ntsr1.evoked_spec=grabSpecInTimeWindow(theta_noLED.allS.S{j},theta_noLED.allS.t,evokedWindow);
                    all_theta_LED_Ntsr1.spont_spec=grabSpecInTimeWindow(theta_LED.allS.S{j},theta_LED.allS.t,spontWindow);
                    all_theta_LED_Ntsr1.evoked_spec=grabSpecInTimeWindow(theta_LED.allS.S{j},theta_LED.allS.t,evokedWindow);
                    
                    [all_noTheta_noLED_Ntsr1.psth,all_noTheta_LED_Ntsr1.psth,all_theta_noLED_Ntsr1.psth,all_theta_LED_Ntsr1.psth,p_noTheta,p_theta]=divideUpPSTH(psth,noThetaTrials,j,spontWindow,evokedWindow);
                    psth_t=psth.t;
                    
                    all_p_noTheta_noLED_Ntsr1=[all_p_noTheta_noLED_Ntsr1 p_noTheta];
                    all_p_theta_noLED_Ntsr1=[all_p_theta_noLED_Ntsr1 p_theta];
                else
                    if isFs(j)==takeFSForNonNtsr1
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
                        
                        all_noTheta_noLED.spont_spec=grabSpecInTimeWindow(noTheta_noLED.allS.S{j},noTheta_noLED.allS.t,spontWindow);
                        all_noTheta_noLED.evoked_spec=grabSpecInTimeWindow(noTheta_noLED.allS.S{j},noTheta_noLED.allS.t,evokedWindow);
                        all_noTheta_LED.spont_spec=grabSpecInTimeWindow(noTheta_LED.allS.S{j},noTheta_LED.allS.t,spontWindow);
                        all_noTheta_LED.evoked_spec=grabSpecInTimeWindow(noTheta_LED.allS.S{j},noTheta_LED.allS.t,evokedWindow);
                        all_theta_noLED.spont_spec=grabSpecInTimeWindow(theta_noLED.allS.S{j},theta_noLED.allS.t,spontWindow);
                        all_theta_noLED.evoked_spec=grabSpecInTimeWindow(theta_noLED.allS.S{j},theta_noLED.allS.t,evokedWindow);
                        all_theta_LED.spont_spec=grabSpecInTimeWindow(theta_LED.allS.S{j},theta_LED.allS.t,spontWindow);
                        all_theta_LED.evoked_spec=grabSpecInTimeWindow(theta_LED.allS.S{j},theta_LED.allS.t,evokedWindow);
                        
                        [all_noTheta_noLED.psth,all_noTheta_LED.psth,all_theta_noLED.psth,all_theta_LED.psth,p_noTheta,p_theta]=divideUpPSTH(psth,noThetaTrials,j,spontWindow,evokedWindow);
                        psth_t=psth.t;
                    end
                end
            else
                if classifyAsNtsr1(j)==1
                    [noTheta_noLED.allS.t,noTheta_noLED.allS.f,noTheta_noLED.allS.S{j}]=cutToSize(all_noTheta_noLED_Ntsr1,noTheta_noLED.allS.t,noTheta_noLED.allS.f,noTheta_noLED.allS.S{j});
                    [noTheta_LED.allS.t,noTheta_LED.allS.f,noTheta_LED.allS.S{j}]=cutToSize(all_noTheta_LED_Ntsr1,noTheta_LED.allS.t,noTheta_LED.allS.f,noTheta_LED.allS.S{j});
                    [theta_noLED.allS.t,theta_noLED.allS.f,theta_noLED.allS.S{j}]=cutToSize(all_theta_noLED_Ntsr1,theta_noLED.allS.t,theta_noLED.allS.f,theta_noLED.allS.S{j});
                    [theta_LED.allS.t,theta_LED.allS.f,theta_LED.allS.S{j}]=cutToSize(all_theta_LED_Ntsr1,theta_LED.allS.t,theta_LED.allS.f,theta_LED.allS.S{j});
                    psth=PSTHcutToSize(all_noTheta_noLED.psth,psth);
                    
                    tmp=cat(3,all_noTheta_noLED_Ntsr1.S,noTheta_noLED.allS.S{j}); 
                    all_noTheta_noLED_Ntsr1.S=nansum(tmp,3);
                    tmp=cat(3,all_noTheta_LED_Ntsr1.S,noTheta_LED.allS.S{j}); 
                    all_noTheta_LED_Ntsr1.S=nansum(tmp,3);
                    tmp=cat(3,all_theta_noLED_Ntsr1.S,theta_noLED.allS.S{j}); 
                    all_theta_noLED_Ntsr1.S=nansum(tmp,3);
                    tmp=cat(3,all_theta_LED_Ntsr1.S,theta_LED.allS.S{j}); 
                    all_theta_LED_Ntsr1.S=nansum(tmp,3);
                    runningTallyNtsr1=runningTallyNtsr1+1;
                    
                    all_noTheta_noLED_Ntsr1.spont_spec=[all_noTheta_noLED_Ntsr1.spont_spec; grabSpecInTimeWindow(noTheta_noLED.allS.S{j},noTheta_noLED.allS.t,spontWindow)];
                    all_noTheta_noLED_Ntsr1.evoked_spec=[all_noTheta_noLED_Ntsr1.evoked_spec; grabSpecInTimeWindow(noTheta_noLED.allS.S{j},noTheta_noLED.allS.t,evokedWindow)];
                    all_noTheta_LED_Ntsr1.spont_spec=[all_noTheta_LED_Ntsr1.spont_spec; grabSpecInTimeWindow(noTheta_LED.allS.S{j},noTheta_LED.allS.t,spontWindow)];
                    all_noTheta_LED_Ntsr1.evoked_spec=[all_noTheta_LED_Ntsr1.evoked_spec; grabSpecInTimeWindow(noTheta_LED.allS.S{j},noTheta_LED.allS.t,evokedWindow)];
                    all_theta_noLED_Ntsr1.spont_spec=[all_theta_noLED_Ntsr1.spont_spec; grabSpecInTimeWindow(theta_noLED.allS.S{j},theta_noLED.allS.t,spontWindow)];
                    all_theta_noLED_Ntsr1.evoked_spec=[all_theta_noLED_Ntsr1.evoked_spec; grabSpecInTimeWindow(theta_noLED.allS.S{j},theta_noLED.allS.t,evokedWindow)];
                    all_theta_LED_Ntsr1.spont_spec=[all_theta_LED_Ntsr1.spont_spec; grabSpecInTimeWindow(theta_LED.allS.S{j},theta_LED.allS.t,spontWindow)];
                    all_theta_LED_Ntsr1.evoked_spec=[all_theta_LED_Ntsr1.evoked_spec; grabSpecInTimeWindow(theta_LED.allS.S{j},theta_LED.allS.t,evokedWindow)];
                    
                    [all_noTheta_noLED_psth,all_noTheta_LED_psth,all_theta_noLED_psth,all_theta_LED_psth,p_noTheta,p_theta]=divideUpPSTH(psth,noThetaTrials,j,spontWindow,evokedWindow);
                    all_noTheta_noLED_Ntsr1.psth=[all_noTheta_noLED_Ntsr1.psth; all_noTheta_noLED_psth];
                    all_noTheta_LED_Ntsr1.psth=[all_noTheta_LED_Ntsr1.psth; all_noTheta_LED_psth];
                    all_theta_noLED_Ntsr1.psth=[all_theta_noLED_Ntsr1.psth; all_theta_noLED_psth];
                    all_theta_LED_Ntsr1.psth=[all_theta_LED_Ntsr1.psth; all_theta_LED_psth];
                    
                    all_p_noTheta_noLED_Ntsr1=[all_p_noTheta_noLED_Ntsr1 p_noTheta];
                    all_p_theta_noLED_Ntsr1=[all_p_theta_noLED_Ntsr1 p_theta];
                    
                    if size(noTheta_noLED.HFa,2)>size(all_noTheta_noLED_Ntsr1.HFa,2)
                        all_noTheta_noLED_Ntsr1.HFa=[all_noTheta_noLED_Ntsr1.HFa; noTheta_noLED.HFa(:,1:size(all_noTheta_noLED_Ntsr1.HFa,2))];
                        all_noTheta_LED_Ntsr1.HFa=[all_noTheta_LED_Ntsr1.HFa; noTheta_LED.HFa(:,1:size(all_noTheta_LED_Ntsr1.HFa,2))];
                        all_theta_noLED_Ntsr1.HFa=[all_theta_noLED_Ntsr1.HFa; theta_noLED.HFa(:,1:size(all_theta_noLED_Ntsr1.HFa,2))];
                        all_theta_LED_Ntsr1.HFa=[all_theta_LED_Ntsr1.HFa; theta_LED.HFa(:,1:size(all_theta_LED_Ntsr1.HFa,2))];
                        
                        all_noTheta_noLED_Ntsr1.LFa=[all_noTheta_noLED_Ntsr1.LFa; noTheta_noLED.LFa(:,1:size(all_noTheta_noLED_Ntsr1.LFa,2))];
                        all_noTheta_LED_Ntsr1.LFa=[all_noTheta_LED_Ntsr1.LFa; noTheta_LED.LFa(:,1:size(all_noTheta_LED_Ntsr1.LFa,2))];
                        all_theta_noLED_Ntsr1.LFa=[all_theta_noLED_Ntsr1.LFa; theta_noLED.LFa(:,1:size(all_theta_noLED_Ntsr1.LFa,2))];
                        all_theta_LED_Ntsr1.LFa=[all_theta_LED_Ntsr1.LFa; theta_LED.LFa(:,1:size(all_theta_LED_Ntsr1.LFa,2))];
                        
                        all_noTheta_noLED_Ntsr1.F1amp=[all_noTheta_noLED_Ntsr1.F1amp; noTheta_noLED.F1amp(:,1:size(all_noTheta_noLED_Ntsr1.F1amp,2))];
                        all_noTheta_LED_Ntsr1.F1amp=[all_noTheta_LED_Ntsr1.F1amp; noTheta_LED.F1amp(:,1:size(all_noTheta_LED_Ntsr1.F1amp,2))];
                        all_theta_noLED_Ntsr1.F1amp=[all_theta_noLED_Ntsr1.F1amp; theta_noLED.F1amp(:,1:size(all_theta_noLED_Ntsr1.F1amp,2))];
                        all_theta_LED_Ntsr1.F1amp=[all_theta_LED_Ntsr1.F1amp; theta_LED.F1amp(:,1:size(all_theta_LED_Ntsr1.F1amp,2))];
                    elseif size(noTheta_noLED.HFa,2)<size(all_noTheta_noLED_Ntsr1.HFa,2)
                        all_noTheta_noLED_Ntsr1.HFa=[all_noTheta_noLED_Ntsr1.HFa; [noTheta_noLED.HFa nan(size(noTheta_noLED.HFa,1),size(all_noTheta_noLED_Ntsr1.HFa,2)-size(noTheta_noLED.HFa,2))]];
                        all_noTheta_LED_Ntsr1.HFa=[all_noTheta_LED_Ntsr1.HFa; [noTheta_LED.HFa nan(size(noTheta_LED.HFa,1),size(all_noTheta_LED_Ntsr1.HFa,2)-size(noTheta_LED.HFa,2))]];
                        all_theta_noLED_Ntsr1.HFa=[all_theta_noLED_Ntsr1.HFa; [theta_noLED.HFa nan(size(theta_noLED.HFa,1),size(all_theta_noLED_Ntsr1.HFa,2)-size(theta_noLED.HFa,2))]];
                        all_theta_LED_Ntsr1.HFa=[all_theta_LED_Ntsr1.HFa; [theta_LED.HFa nan(size(theta_LED.HFa,1),size(all_theta_LED_Ntsr1.HFa,2)-size(theta_LED.HFa,2))]];
                        
                        all_noTheta_noLED_Ntsr1.LFa=[all_noTheta_noLED_Ntsr1.LFa; [noTheta_noLED.LFa nan(size(noTheta_noLED.LFa,1),size(all_noTheta_noLED_Ntsr1.LFa,2)-size(noTheta_noLED.LFa,2))]];
                        all_noTheta_LED_Ntsr1.LFa=[all_noTheta_LED_Ntsr1.LFa; [noTheta_LED.LFa nan(size(noTheta_LED.LFa,1),size(all_noTheta_LED_Ntsr1.LFa,2)-size(noTheta_LED.LFa,2))]];
                        all_theta_noLED_Ntsr1.LFa=[all_theta_noLED_Ntsr1.LFa; [theta_noLED.LFa nan(size(theta_noLED.LFa,1),size(all_theta_noLED_Ntsr1.LFa,2)-size(theta_noLED.LFa,2))]];
                        all_theta_LED_Ntsr1.LFa=[all_theta_LED_Ntsr1.LFa; [theta_LED.LFa nan(size(theta_LED.LFa,1),size(all_theta_LED_Ntsr1.LFa,2)-size(theta_LED.LFa,2))]];
                        
                        all_noTheta_noLED_Ntsr1.F1amp=[all_noTheta_noLED_Ntsr1.F1amp; [noTheta_noLED.F1amp nan(size(noTheta_noLED.F1amp,1),size(all_noTheta_noLED_Ntsr1.F1amp,2)-size(noTheta_noLED.F1amp,2))]];
                        all_noTheta_LED_Ntsr1.F1amp=[all_noTheta_LED_Ntsr1.F1amp; [noTheta_LED.F1amp nan(size(noTheta_LED.F1amp,1),size(all_noTheta_LED_Ntsr1.F1amp,2)-size(noTheta_LED.F1amp,2))]];
                        all_theta_noLED_Ntsr1.F1amp=[all_theta_noLED_Ntsr1.F1amp; [theta_noLED.F1amp nan(size(theta_noLED.F1amp,1),size(all_theta_noLED_Ntsr1.F1amp,2)-size(theta_noLED.F1amp,2))]];
                        all_theta_LED_Ntsr1.F1amp=[all_theta_LED_Ntsr1.F1amp; [theta_LED.F1amp nan(size(theta_LED.F1amp,1),size(all_theta_LED_Ntsr1.F1amp,2)-size(theta_LED.F1amp,2))]];
                    end   
                else
                    if isFs(j)==takeFSForNonNtsr1
                        [noTheta_noLED.allS.t,noTheta_noLED.allS.f,noTheta_noLED.allS.S{j}]=cutToSize(all_noTheta_noLED,noTheta_noLED.allS.t,noTheta_noLED.allS.f,noTheta_noLED.allS.S{j});
                        [noTheta_LED.allS.t,noTheta_LED.allS.f,noTheta_LED.allS.S{j}]=cutToSize(all_noTheta_LED,noTheta_LED.allS.t,noTheta_LED.allS.f,noTheta_LED.allS.S{j});
                        [theta_noLED.allS.t,theta_noLED.allS.f,theta_noLED.allS.S{j}]=cutToSize(all_theta_noLED,theta_noLED.allS.t,theta_noLED.allS.f,theta_noLED.allS.S{j});
                        [theta_LED.allS.t,theta_LED.allS.f,theta_LED.allS.S{j}]=cutToSize(all_theta_LED,theta_LED.allS.t,theta_LED.allS.f,theta_LED.allS.S{j});
                        psth=PSTHcutToSize(all_noTheta_noLED.psth,psth);
                        
                        tmp=cat(3,all_noTheta_noLED.S,noTheta_noLED.allS.S{j});
                        all_noTheta_noLED.S=nansum(tmp,3);
                        tmp=cat(3,all_noTheta_LED.S,noTheta_LED.allS.S{j});
                        all_noTheta_LED.S=nansum(tmp,3);
                        tmp=cat(3,all_theta_noLED.S,theta_noLED.allS.S{j});
                        all_theta_noLED.S=nansum(tmp,3);
                        tmp=cat(3,all_theta_LED.S,theta_LED.allS.S{j});
                        all_theta_LED.S=nansum(tmp,3);
                        runningTallyOther=runningTallyOther+1;
                        
                        all_noTheta_noLED.spont_spec=[all_noTheta_noLED.spont_spec; grabSpecInTimeWindow(noTheta_noLED.allS.S{j},noTheta_noLED.allS.t,spontWindow)];
                        all_noTheta_noLED.evoked_spec=[all_noTheta_noLED.evoked_spec; grabSpecInTimeWindow(noTheta_noLED.allS.S{j},noTheta_noLED.allS.t,evokedWindow)];
                        all_noTheta_LED.spont_spec=[all_noTheta_LED.spont_spec; grabSpecInTimeWindow(noTheta_LED.allS.S{j},noTheta_LED.allS.t,spontWindow)];
                        all_noTheta_LED.evoked_spec=[all_noTheta_LED.evoked_spec; grabSpecInTimeWindow(noTheta_LED.allS.S{j},noTheta_LED.allS.t,evokedWindow)];
                        all_theta_noLED.spont_spec=[all_theta_noLED.spont_spec; grabSpecInTimeWindow(theta_noLED.allS.S{j},theta_noLED.allS.t,spontWindow)];
                        all_theta_noLED.evoked_spec=[all_theta_noLED.evoked_spec; grabSpecInTimeWindow(theta_noLED.allS.S{j},theta_noLED.allS.t,evokedWindow)];
                        all_theta_LED.spont_spec=[all_theta_LED.spont_spec; grabSpecInTimeWindow(theta_LED.allS.S{j},theta_LED.allS.t,spontWindow)];
                        all_theta_LED.evoked_spec=[all_theta_LED.evoked_spec; grabSpecInTimeWindow(theta_LED.allS.S{j},theta_LED.allS.t,evokedWindow)];
                        
                        [all_noTheta_noLED_psth,all_noTheta_LED_psth,all_theta_noLED_psth,all_theta_LED_psth,p_noTheta,p_theta]=divideUpPSTH(psth,noThetaTrials,j,spontWindow,evokedWindow);
                        all_noTheta_noLED.psth=[all_noTheta_noLED.psth; all_noTheta_noLED_psth];
                        all_noTheta_LED.psth=[all_noTheta_LED.psth; all_noTheta_LED_psth];
                        all_theta_noLED.psth=[all_theta_noLED.psth; all_theta_noLED_psth];
                        all_theta_LED.psth=[all_theta_LED.psth; all_theta_LED_psth];
                        
                        if size(noTheta_noLED.HFa,2)>size(all_noTheta_noLED.HFa,2)
                            all_noTheta_noLED.HFa=[all_noTheta_noLED.HFa; noTheta_noLED.HFa(:,1:size(all_noTheta_noLED.HFa,2))];
                            all_noTheta_LED.HFa=[all_noTheta_LED.HFa; noTheta_LED.HFa(:,1:size(all_noTheta_LED.HFa,2))];
                            all_theta_noLED.HFa=[all_theta_noLED.HFa; theta_noLED.HFa(:,1:size(all_theta_noLED.HFa,2))];
                            all_theta_LED.HFa=[all_theta_LED.HFa; theta_LED.HFa(:,1:size(all_theta_LED.HFa,2))];
                            
                            all_noTheta_noLED.LFa=[all_noTheta_noLED.LFa; noTheta_noLED.LFa(:,1:size(all_noTheta_noLED.LFa,2))];
                            all_noTheta_LED.LFa=[all_noTheta_LED.LFa; noTheta_LED.LFa(:,1:size(all_noTheta_LED.LFa,2))];
                            all_theta_noLED.LFa=[all_theta_noLED.LFa; theta_noLED.LFa(:,1:size(all_theta_noLED.LFa,2))];
                            all_theta_LED.LFa=[all_theta_LED.LFa; theta_LED.LFa(:,1:size(all_theta_LED.LFa,2))];

                            all_noTheta_noLED.F1amp=[all_noTheta_noLED.F1amp; noTheta_noLED.F1amp(:,1:size(all_noTheta_noLED.F1amp,2))];
                            all_noTheta_LED.F1amp=[all_noTheta_LED.F1amp; noTheta_LED.F1amp(:,1:size(all_noTheta_LED.F1amp,2))];
                            all_theta_noLED.F1amp=[all_theta_noLED.F1amp; theta_noLED.F1amp(:,1:size(all_theta_noLED.F1amp,2))];
                            all_theta_LED.F1amp=[all_theta_LED.F1amp; theta_LED.F1amp(:,1:size(all_theta_LED.F1amp,2))];
                        elseif size(noTheta_noLED.HFa,2)<size(all_noTheta_noLED.HFa,2)
                            all_noTheta_noLED.HFa=[all_noTheta_noLED.HFa; [noTheta_noLED.HFa nan(size(noTheta_noLED.HFa,1),size(all_noTheta_noLED.HFa,2)-size(noTheta_noLED.HFa,2))]];
                            all_noTheta_LED.HFa=[all_noTheta_LED.HFa; [noTheta_LED.HFa nan(size(noTheta_LED.HFa,1),size(all_noTheta_LED.HFa,2)-size(noTheta_LED.HFa,2))]];
                            all_theta_noLED.HFa=[all_theta_noLED.HFa; [theta_noLED.HFa nan(size(theta_noLED.HFa,1),size(all_theta_noLED.HFa,2)-size(theta_noLED.HFa,2))]];
                            all_theta_LED.HFa=[all_theta_LED.HFa; [theta_LED.HFa nan(size(theta_LED.HFa,1),size(all_theta_LED.HFa,2)-size(theta_LED.HFa,2))]];
                            
                            all_noTheta_noLED.LFa=[all_noTheta_noLED.LFa; [noTheta_noLED.LFa nan(size(noTheta_noLED.LFa,1),size(all_noTheta_noLED.LFa,2)-size(noTheta_noLED.LFa,2))]];
                            all_noTheta_LED.LFa=[all_noTheta_LED.LFa; [noTheta_LED.LFa nan(size(noTheta_LED.LFa,1),size(all_noTheta_LED.LFa,2)-size(noTheta_LED.LFa,2))]];
                            all_theta_noLED.LFa=[all_theta_noLED.LFa; [theta_noLED.LFa nan(size(theta_noLED.LFa,1),size(all_theta_noLED.LFa,2)-size(theta_noLED.LFa,2))]];
                            all_theta_LED.LFa=[all_theta_LED.LFa; [theta_LED.LFa nan(size(theta_LED.LFa,1),size(all_theta_LED.LFa,2)-size(theta_LED.LFa,2))]];
                            
                            all_noTheta_noLED.F1amp=[all_noTheta_noLED.F1amp; [noTheta_noLED.F1amp nan(size(noTheta_noLED.F1amp,1),size(all_noTheta_noLED.F1amp,2)-size(noTheta_noLED.F1amp,2))]];
                            all_noTheta_LED.F1amp=[all_noTheta_LED.F1amp; [noTheta_LED.F1amp nan(size(noTheta_LED.F1amp,1),size(all_noTheta_LED.F1amp,2)-size(noTheta_LED.F1amp,2))]];
                            all_theta_noLED.F1amp=[all_theta_noLED.F1amp; [theta_noLED.F1amp nan(size(theta_noLED.F1amp,1),size(all_theta_noLED.F1amp,2)-size(theta_noLED.F1amp,2))]];
                            all_theta_LED.F1amp=[all_theta_LED.F1amp; [theta_LED.F1amp nan(size(theta_LED.F1amp,1),size(all_theta_LED.F1amp,2)-size(theta_LED.F1amp,2))]];
                        end
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

ds=1;
y1_noTheta=plotWStderr(downSampAv(all_noTheta_noLED_Ntsr1.t,ds),downSampMatrix(all_noTheta_noLED_Ntsr1.F1amp,ds),[],'k','c',1);
title('No theta Ntsr1 F1 amp');

ds=1;
y1_theta=plotWStderr(downSampAv(all_noTheta_noLED_Ntsr1.t,ds),downSampMatrix(all_theta_noLED_Ntsr1.F1amp,ds),[],'k','c',1);
title('Theta Ntsr1 F1 amp');

ds=5;
plotWStderr(downSampAv(psth_t,ds),downSampMatrix(all_noTheta_noLED_Ntsr1.psth,ds),[],'k','c',1);
title('No theta Ntsr1');

plotWStderr(downSampAv(psth_t,ds),downSampMatrix(all_theta_noLED_Ntsr1.psth,ds),[],'k','c',1);
title('Theta Ntsr1');

ev_noTheta=nanmean(all_noTheta_noLED_Ntsr1.psth(:,psth_t>=evokedWindow(1) & psth_t<=evokedWindow(2)),2);
spont_noTheta=nanmean(all_noTheta_noLED_Ntsr1.psth(:,psth_t>=spontWindow(1) & psth_t<=spontWindow(2)),2);
figure();
for i=1:length(spont_noTheta)
    if all_p_noTheta_noLED_Ntsr1(i)<0.05
        line([1 2],[spont_noTheta(i) ev_noTheta(i)],'Color','k','LineWidth',1.8);
    else
        line([1 2],[spont_noTheta(i) ev_noTheta(i)],'Color','k','LineWidth',0.5);
    end
    hold on;
end
set(gca,'Yscale','log');
disp('p-value for raw firing rates spont vs evoked no theta');
disp(signrank(ev_noTheta,spont_noTheta));

ev_theta=nanmean(all_theta_noLED_Ntsr1.psth(:,psth_t>=evokedWindow(1) & psth_t<=evokedWindow(2)),2);
spont_theta=nanmean(all_theta_noLED_Ntsr1.psth(:,psth_t>=spontWindow(1) & psth_t<=spontWindow(2)),2);
figure();
for i=1:length(spont_theta)
    if all_p_theta_noLED_Ntsr1(i)<0.05
        line([1 2],[spont_theta(i) ev_theta(i)],'Color','k','LineWidth',1.8);
    else
        line([1 2],[spont_theta(i) ev_theta(i)],'Color','k','LineWidth',0.5);
    end
    hold on;
end
set(gca,'Yscale','log');
disp('p-value for raw firing rates spont vs evoked theta');
disp(signrank(ev_theta,spont_theta));

figure(); 
scatter(ev_noTheta-spont_noTheta,ev_theta-spont_theta,[],'k');
hold on;
scatter(ev_noTheta(all_p_noTheta_noLED_Ntsr1<0.05 & all_p_theta_noLED_Ntsr1<0.05)-spont_noTheta(all_p_noTheta_noLED_Ntsr1<0.05 & all_p_theta_noLED_Ntsr1<0.05),ev_theta(all_p_noTheta_noLED_Ntsr1<0.05 & all_p_theta_noLED_Ntsr1<0.05)-spont_theta(all_p_noTheta_noLED_Ntsr1<0.05 & all_p_theta_noLED_Ntsr1<0.05),[],'k','filled');
scatter(ev_noTheta(all_p_noTheta_noLED_Ntsr1<0.05 & all_p_theta_noLED_Ntsr1>=0.05)-spont_noTheta(all_p_noTheta_noLED_Ntsr1<0.05 & all_p_theta_noLED_Ntsr1>=0.05),ev_theta(all_p_noTheta_noLED_Ntsr1<0.05 & all_p_theta_noLED_Ntsr1>=0.05)-spont_theta(all_p_noTheta_noLED_Ntsr1<0.05 & all_p_theta_noLED_Ntsr1>=0.05),[],[0.7765 0.2275 0.8863],'filled');
scatter(ev_noTheta(all_p_noTheta_noLED_Ntsr1>=0.05 & all_p_theta_noLED_Ntsr1<0.05)-spont_noTheta(all_p_noTheta_noLED_Ntsr1>=0.05 & all_p_theta_noLED_Ntsr1<0.05),ev_theta(all_p_noTheta_noLED_Ntsr1>=0.05 & all_p_theta_noLED_Ntsr1<0.05)-spont_theta(all_p_noTheta_noLED_Ntsr1>=0.05 & all_p_theta_noLED_Ntsr1<0.05),[],[0 0.7529 0.5529],'filled');

temp1=nanmean(all_noTheta_noLED_Ntsr1.spont_spec(:,all_noTheta_noLED_Ntsr1.f>=9.5 & all_noTheta_noLED_Ntsr1.f<=16.5),2)./nanmean(all_noTheta_noLED_Ntsr1.spont_spec(:,all_noTheta_noLED_Ntsr1.f>=4 & all_noTheta_noLED_Ntsr1.f<=6.5),2);
temp2=nanmean(all_noTheta_noLED_Ntsr1.evoked_spec(:,all_noTheta_noLED_Ntsr1.f>=9.5 & all_noTheta_noLED_Ntsr1.f<=16.5),2)./nanmean(all_noTheta_noLED_Ntsr1.evoked_spec(:,all_noTheta_noLED_Ntsr1.f>=4 & all_noTheta_noLED_Ntsr1.f<=6.5),2);
plotComparisonLines(temp1,temp2);

temp1=nanmean(all_theta_noLED_Ntsr1.spont_spec(:,all_theta_noLED_Ntsr1.f>=9.5 & all_theta_noLED_Ntsr1.f<=16.5),2)./nanmean(all_theta_noLED_Ntsr1.spont_spec(:,all_theta_noLED_Ntsr1.f>=4 & all_theta_noLED_Ntsr1.f<=6.5),2);
temp2=nanmean(all_theta_noLED_Ntsr1.evoked_spec(:,all_theta_noLED_Ntsr1.f>=9.5 & all_theta_noLED_Ntsr1.f<=16.5),2)./nanmean(all_theta_noLED_Ntsr1.evoked_spec(:,all_theta_noLED_Ntsr1.f>=4 & all_theta_noLED_Ntsr1.f<=6.5),2);
plotComparisonLines(temp1,temp2);

temp1=nanmean(all_noTheta_noLED.spont_spec(:,all_noTheta_noLED.f>=9.5 & all_noTheta_noLED.f<=16.5),2)./nanmean(all_noTheta_noLED.spont_spec(:,all_noTheta_noLED.f>=4 & all_noTheta_noLED.f<=6.5),2);
temp2=nanmean(all_noTheta_noLED.evoked_spec(:,all_noTheta_noLED.f>=9.5 & all_noTheta_noLED.f<=16.5),2)./nanmean(all_noTheta_noLED.evoked_spec(:,all_noTheta_noLED.f>=4 & all_noTheta_noLED.f<=6.5),2);
plotComparisonLines(temp1,temp2);

% temp=all_noTheta_noLED_Ntsr1.spont_spec(:,all_noTheta_noLED_Ntsr1.f>=4.688 & all_noTheta_noLED_Ntsr1.f<=50);
temp=all_noTheta_noLED_Ntsr1.spont_spec(:,all_noTheta_noLED_Ntsr1.f>=0 & all_noTheta_noLED_Ntsr1.f<=50);
temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
all_noTheta_noLED_Ntsr1.spont_spec=temp;
temp=all_noTheta_noLED_Ntsr1.evoked_spec(:,all_noTheta_noLED_Ntsr1.f>=0 & all_noTheta_noLED_Ntsr1.f<=50);
temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
all_noTheta_noLED_Ntsr1.evoked_spec=temp;
plotWStderr(all_noTheta_noLED_Ntsr1.f(all_noTheta_noLED_Ntsr1.f>=0 & all_noTheta_noLED_Ntsr1.f<=50),all_noTheta_noLED_Ntsr1.spont_spec,all_noTheta_noLED_Ntsr1.evoked_spec,'k','r',1);
title('Ntsr1 No theta No LED spont vs evoked');
set(gca,'xscale','log');

temp=all_noTheta_noLED.spont_spec(:,all_noTheta_noLED.f>=0 & all_noTheta_noLED.f<=50);
temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
all_noTheta_noLED.spont_spec=temp;
temp=all_noTheta_noLED.evoked_spec(:,all_noTheta_noLED.f>=0 & all_noTheta_noLED.f<=50);
temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
all_noTheta_noLED.evoked_spec=temp;
plotWStderr(all_noTheta_noLED.f(all_noTheta_noLED.f>=0 & all_noTheta_noLED.f<=50),all_noTheta_noLED.spont_spec,all_noTheta_noLED.evoked_spec,'k','r',1);
title('No theta No LED spont vs evoked');
set(gca,'xscale','log');



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
temp=all_noTheta_noLED.S(~isnan(all_noTheta_noLED.t) & all_noTheta_noLED.t<=cutToMaxTime & all_noTheta_noLED.t>=cutToMinTime,all_noTheta_noLED.f<=50);
temp=temp./repmat(nansum(temp,2),1,size(temp,2));
imagesc(all_noTheta_noLED.t(~isnan(all_noTheta_noLED.t) & all_noTheta_noLED.t<=cutToMaxTime & all_noTheta_noLED.t>=cutToMinTime),all_noTheta_noLED.f(all_noTheta_noLED.f<=50),temp');
title('noTheta noLED normalized');

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
current_t=all_noTheta_noLED_Ntsr1.t(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime);
temp=all_noTheta_noLED_Ntsr1.S(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime,all_noTheta_noLED_Ntsr1.f<=50);
temp(current_t>3 & current_t<=3.725,:)=flipud(0.7*temp(current_t>2.3 & current_t<=3.025,:))+0.3*temp(current_t>3.3 & current_t<=4.025,:);
temp(:,1:3)=temp(:,1:3).*0.17.*(0.7225/0.1121);
temp(:,4)=temp(:,4).*0.73.*(0.7225/0.5035);
temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
imagesc(downSampAv(all_noTheta_noLED_Ntsr1.t(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),1),all_noTheta_noLED_Ntsr1.f(all_noTheta_noLED_Ntsr1.f<=50),downSampMatrix(temp',1));
title('noTheta noLED Ntsr1 normalized');

dataInt=interp2(temp,4);
figure();
imagesc(all_noTheta_noLED_Ntsr1.t(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime),all_noTheta_noLED_Ntsr1.f(all_noTheta_noLED_Ntsr1.f<=50),dataInt');
pause();

% tempie=[nan(4,size(allVisBaseSubByFreq_noTheta,2)+8); [nan(size(allVisBaseSubByFreq_noTheta,1),4) allVisBaseSubByFreq_noTheta nan(size(allVisBaseSubByFreq_noTheta,1),4)]; nan(4,size(allVisBaseSubByFreq_noTheta,2)+8)];
% K=ones([5 5]); figure(); imagesc(interp2(conv2(tempie,K),4)');
% temp=conv2(tempie,K);
% temp=temp(1:end-1,1:end-1);
% temp=interp2(temp,4);
% 
% F1power_all=temp;
% 
% tempie=[nan(4,size(p_noTheta',2)+8); [nan(size(p_noTheta',1),4) p_noTheta' nan(size(p_noTheta',1),4)]; nan(4,size(p_noTheta',2)+8)];
% K=ones([5 5]); figure(); imagesc(interp2(conv2(tempie,K),4)');
% temp=conv2(tempie,K);
% temp=temp(1:end-1,1:end-1);
% temp=interp2(temp,4);
% 
% p_all=temp;
% 
% temperfid=F1power_all; 
% % temperfid(p_all>3*4.1532)=nanmean(nanmean(F1power_all)); % no theta
% % temperfid(p_all>23.1375)=nanmean(nanmean(F1power_all)); % theta
% temperfid(p_all>15.16)=nanmean(nanmean(F1power_all)); % no theta
% figure(); 
% imagesc(conv2(temperfid,K,'same')');
% 
% temp=temperfid;
% 
% bin1_temp=nanmean(temp(112+1:112+12,:),1);
% bin2_temp=nanmean(temp(112+13:112+13+12,:),1);
% bin3_temp=nanmean(temp(112+13+12:112+13+2*12,:),1);
% bin4_temp=nanmean(temp(112+13+2*12:112+13+3*12,:),1);
% bin5_temp=nanmean(temp(112+13+3*12:112+13+4*12,:),1);
% bin6_temp=nanmean(temp(112+13+4*12:112+13+5*12,:),1);
% bin7_temp=nanmean(temp(112+13+5*12:112+13+6*12,:),1);
% bin8_temp=nanmean(temp(112+13+6*12:112+13+7*12,:),1);
% bin9_temp=nanmean(temp(112+13+7*12:112+13+8*12,:),1);
% bin10_temp=nanmean(temp(112+13+8*12:112+13+9*12,:),1);
% bin11_temp=nanmean(temp(112+13+9*12:112+13+10*12,:),1);
% bin12_temp=nanmean(temp(112+13+10*12:112+13+11*12,:),1);
% bin13_temp=nanmean(temp(112+13+11*12:112+13+12*12,:),1);
% bin14_temp=nanmean(temp(112+13+12*12:112+13+13*12,:),1);
% bin15_temp=nanmean(temp(112+13+13*12:end-112-1,:),1);
% 
% all(1:12,:)=repmat(bin1_temp,12,1);
% all(13:13+12-1,:)=repmat(bin2_temp,12,1);
% all(13+12:13+12+8-1,:)=repmat(bin3_temp,8,1);
% all(13+12+8:13+12+8+6-1,:)=repmat(bin4_temp,6,1);
% all(13+12+8+6:13+12+8+6+2-1,:)=repmat(bin5_temp,2,1);
% all(13+12+8+6+2:13+12+8+6+2+4-1,:)=repmat(bin6_temp,4,1);
% all(13+12+8+6+2+4:13+12+8+6+2+4+3-1,:)=repmat(bin7_temp,3,1);
% all(13+12+8+6+2+4+3:13+12+8+6+2+4+3+2-1,:)=repmat(bin8_temp,2,1);
% all(13+12+8+6+2+4+3+2:13+12+8+6+2+4+3+2+1-1,:)=repmat(bin9_temp,1,1);
% all(13+12+8+6+2+4+3+2+1:13+12+8+6+2+4+3+2+1+1-1,:)=repmat(bin10_temp,1,1);
% all(13+12+8+6+2+4+3+2+1+1:13+12+8+6+2+4+3+2+1+1+7-1,:)=repmat(bin11_temp,7,1);
% all(13+12+8+6+2+4+3+2+1+1+7:13+12+8+6+2+4+3+2+1+1+7+5-1,:)=repmat(bin12_temp,5,1);
% all(13+12+8+6+2+4+3+2+1+1+7+5:13+12+8+6+2+4+3+2+1+1+7+5+4-1,:)=repmat(bin13_temp,4,1);
% all(13+12+8+6+2+4+3+2+1+1+7+5+4:13+12+8+6+2+4+3+2+1+1+7+5+4+3-1,:)=repmat(bin14_temp,3,1);
% all(13+12+8+6+2+4+3+2+1+1+7+5+4+3:13+12+8+6+2+4+3+2+1+1+7+5+4+3+2-1,:)=repmat(bin15_temp,2,1);



% all(1,:)=bin1_temp;
% all(2,:)=bin2_temp;
% all(3,:)=bin2_temp;
% all(4,:)=bin3_temp;
% all(5,:)=bin3_temp;
% all(6,:)=bin4_temp;
% all(7,:)=bin4_temp;
% all(8,:)=bin5_temp;
% all(9,:)=bin5_temp;
% all(10,:)=bin6_temp;
% all(11,:)=bin6_temp;
% all(12,:)=bin7_temp;
% all(13,:)=bin7_temp;
% all(14,:)=bin8_temp;
% all(15,:)=bin8_temp;
% all(16,:)=bin9_temp;
% all(17,:)=bin9_temp;
% all(18,:)=bin10_temp;
% all(19,:)=bin10_temp;
% all(20,:)=bin11_temp;
% all(21,:)=bin11_temp;
% all(22,:)=bin11_temp;
% all(23,:)=bin11_temp;
% all(24,:)=bin11_temp;
% all(25,:)=bin11_temp;
% all(26,:)=bin11_temp;
% all(27,:)=bin11_temp;
% all(28,:)=bin11_temp;
% all(29,:)=bin11_temp;
% all(30,:)=bin12_temp;
% all(31,:)=bin12_temp;
% all(32,:)=bin12_temp;
% all(33,:)=bin12_temp;
% all(34,:)=bin12_temp;
% all(35,:)=bin12_temp;
% all(36,:)=bin12_temp;
% all(37,:)=bin12_temp;
% all(38,:)=bin12_temp;
% all(39,:)=bin12_temp;
% all(40,:)=bin12_temp;

figure(); 
temp=all_theta_noLED_Ntsr1.S(~isnan(all_theta_noLED_Ntsr1.t) & all_theta_noLED_Ntsr1.t<=cutToMaxTime & all_theta_noLED_Ntsr1.t>=cutToMinTime,all_theta_noLED_Ntsr1.f<=50);
temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
imagesc(downSampAv(all_theta_noLED_Ntsr1.t(~isnan(all_theta_noLED_Ntsr1.t) & all_theta_noLED_Ntsr1.t<=cutToMaxTime & all_theta_noLED_Ntsr1.t>=cutToMinTime),1),all_theta_noLED_Ntsr1.f(all_theta_noLED_Ntsr1.f<=50),downSampMatrix(temp',1));
title('theta noLED Ntsr1 normalized');

figure(); 
temp=all_noTheta_noLED.S(~isnan(all_noTheta_noLED.t) & all_noTheta_noLED.t<=cutToMaxTime & all_noTheta_noLED.t>=cutToMinTime,all_noTheta_noLED.f<=50);
temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
imagesc(all_noTheta_noLED.t(~isnan(all_noTheta_noLED.t) & all_noTheta_noLED.t<=cutToMaxTime & all_noTheta_noLED.t>=cutToMinTime),all_noTheta_noLED.f(all_noTheta_noLED.f<=50),temp');
title('noTheta noLED normalized');

% x=all_noTheta_noLED_Ntsr1.t(~isnan(all_noTheta_noLED_Ntsr1.t) & all_noTheta_noLED_Ntsr1.t<=cutToMaxTime & all_noTheta_noLED_Ntsr1.t>=cutToMinTime);
% y=all_noTheta_noLED_Ntsr1.f(all_noTheta_noLED_Ntsr1.f<=50);
% smoothBy=3; 
% K=ones(smoothBy);
% smoothMat=conv2(temp,K,'same'); 
% smoothMat=smoothMat(smoothBy+1:end-smoothBy-1,smoothBy+1:end-smoothBy-1);
% figure();
% imagesc(x,y,smoothMat');

end

function plotComparisonLines(grp1,grp2)

figure();
for i=1:length(grp1)
    line([1 2],[grp1(i) grp2(i)]);
    hold on;
end
set(gca,'yscale','log');
p=signrank(grp1,grp2);
disp('p-value');
disp(p);

end

function spec=grabSpecInTimeWindow(specgram,specgram_t,timeWindow)

spec=nanmean(specgram(specgram_t>=timeWindow(1) & specgram_t<=timeWindow(2),:),1);

end

function [y1,y2]=plotWStderr(x,y1,y2,c1,c2,doFill)

baseSub=false;
baseWindow=[1.485 3.535];

if baseSub==true
    y1=y1-repmat(nanmean(y1(:,x>=baseWindow(1) & x<=baseWindow(2)),2),1,size(y1,2));
    if ~isempty(y2)
        y2=y2-repmat(nanmean(y2(:,x>=baseWindow(1) & x<=baseWindow(2)),2),1,size(y2,2));
    end
end

figure();
plot(x,nanmean(y1,1),'Color',c1);
hold on;
if doFill==1
    fill([x fliplr(x)],[nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)) fliplr(nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)))],[0.5 0.5 0.5]);
end
plot(x,nanmean(y1,1),'Color',c1);
plot(x,nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);
plot(x,nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);

if ~isempty(y2)
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
    
end

function [noTheta_noLED,noTheta_LED,theta_noLED,theta_LED,p_noTheta_noLED,p_theta_noLED]=divideUpPSTH(psth,noThetaTrials,takeUnitN,spontWindow,evokedWindow)

l=psth.unitLED{1};
ledOff=0;
ledOn=unique(l(l>0.01));
tri=psth.unitTrials{1};

noTheta_noLED=nan(length(psth.psths),size(psth.psths{1},2));
noTheta_LED=nan(length(psth.psths),size(psth.psths{1},2));
theta_noLED=nan(length(psth.psths),size(psth.psths{1},2));
theta_LED=nan(length(psth.psths),size(psth.psths{1},2));
for i=1:length(psth.psths)
    temp=psth.psths{i};
    noTheta_noLED(i,:)=nanmean(temp(noThetaTrials'==1 & ismember(l,ledOff),:),1);
    psth_t=linspace(nanmin(psth.t),nanmax(psth.t),size(psth.psths{1},2));
    tempie1=nanmean(temp(noThetaTrials'==1 & ismember(l,ledOff),psth_t>=spontWindow(1) & psth_t<=spontWindow(2)),2);
    tempie2=nanmean(temp(noThetaTrials'==1 & ismember(l,ledOff),psth_t>=evokedWindow(1) & psth_t<=evokedWindow(2)),2);
    if all(isnan(tempie1) | isnan(tempie2))
        p_noTheta_noLED(i)=nan;
    else
        p_noTheta_noLED(i)=signrank(tempie1,tempie2);
    end
    noTheta_LED(i,:)=nanmean(temp(noThetaTrials'==1 & ismember(l,ledOn),:),1);
    theta_noLED(i,:)=nanmean(temp(noThetaTrials'==0 & ismember(l,ledOff),:),1);
    tempie1=nanmean(temp(noThetaTrials'==0 & ismember(l,ledOff),psth_t>=spontWindow(1) & psth_t<=spontWindow(2)),2);
    tempie2=nanmean(temp(noThetaTrials'==0 & ismember(l,ledOff),psth_t>=evokedWindow(1) & psth_t<=evokedWindow(2)),2);
    if all(isnan(tempie1) | isnan(tempie2))
        p_theta_noLED(i)=nan;
    else
        p_theta_noLED(i)=signrank(tempie1,tempie2);
    end
    theta_LED(i,:)=nanmean(temp(noThetaTrials'==0 & ismember(l,ledOn),:),1);
end

noTheta_noLED=noTheta_noLED(takeUnitN,:);
noTheta_LED=noTheta_LED(takeUnitN,:);
theta_noLED=theta_noLED(takeUnitN,:);
theta_LED=theta_LED(takeUnitN,:);
p_noTheta_noLED=p_noTheta_noLED(takeUnitN);
p_theta_noLED=p_theta_noLED(takeUnitN);

end

function psthAll=PSTHcutToSize(fitTo,psthAll)

for i=1:length(psthAll.psths)
    psth=psthAll.psths{i};
    if size(psth,2)>size(fitTo,2)
        psth=psth(:,1:size(fitTo,2));
    elseif size(psth,2)<size(fitTo,2)
        psth=[psth nan(size(psth,1),size(fitTo,2)-size(psth,2))];
    end
    psthAll.psths{i}=psth;
end

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

function [noTheta_noLED,noTheta_LED,theta_noLED,theta_LED,psth]=alignToStimOnset(noTheta_noLED,noTheta_LED,theta_noLED,theta_LED,stimWindow,placeStimOnsetAt,psth)

stimOnset=stimWindow(1);
if abs(stimOnset-placeStimOnsetAt)<0.01
    return
end

noTheta_noLED=alignEach(noTheta_noLED,stimOnset,placeStimOnsetAt);
noTheta_LED=alignEach(noTheta_LED,stimOnset,placeStimOnsetAt);
theta_noLED=alignEach(theta_noLED,stimOnset,placeStimOnsetAt);
theta_LED=alignEach(theta_LED,stimOnset,placeStimOnsetAt);

for i=1:length(psth.psths)
    psth.psths{i}=alignPSTH(psth.psths{i},stimOnset,placeStimOnsetAt,psth.t);
end

end

function psth=alignPSTH(psth,stimOnset,placeStimOnsetAt,t)

% assume same binning for all
shiftByInds=floor((stimOnset-placeStimOnsetAt)/(t(2)-t(1)));

if shiftByInds>0
    % shift forward by shiftByInds
    psth=[psth(:,shiftByInds:end) nan(size(psth,1),shiftByInds-1)];
elseif shiftByInds<0
    % shift backward by shiftByInds
    shiftByInds=-shiftByInds;
    psth=[nan(size(psth,1),shiftByInds) psth(:,1:end-shiftByInds)];
end

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
    shiftByInds=-shiftByInds;
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