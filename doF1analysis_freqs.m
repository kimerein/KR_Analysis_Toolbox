function doF1analysis_freqs(expt,useFileInd,outputDir,LFP_Fs,dLGNpsth,usel_noLED,usel_LED,LFPdata,uses,uses_tri,noThetaTrials)

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% Mawake328
% a=1:656;
% uses=   {{[1 1.05]}; {[2 2.05]}; {[4 4.05]}; {[6 6.05]}; {[8 8.05]}; {[10 10.05]}; {[12 12.05]}; {[14 14.05]}; {[16 16.05]}; {[18 18.05]}; {[20 20.05]}; {[30 30.05]}; {[40 40.05]}; {[50 50.05]}; {[60 60.05]}};           
% uses_tri={{a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}};
 
% Mawake327
a=1:913;
uses=   {{[1 1.05]}; {[2 2.05]}; {[4 4.05]}; {[6 6.05]}; {[8 8.05]}; {[10 10.05]}; {[12 12.05]}; {[14 14.05]}; {[16 16.05]}; {[18 18.05]}; {[20 20.05]}; {[30 30.05]}; {[40 40.05]}; {[50 50.05]}; {[60 60.05]}};           
uses_tri={{a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}};

% Mawake329
% a=1:764;
% uses=   {{[1 1.05]}; {[2 2.05]}; {[4 4.05]}; {[6 6.05]}; {[8 8.05]}; {[10 10.05]}; {[12 12.05]}; {[14 14.05]}; {[16 16.05]}; {[18 18.05]}; {[20 20.05]}; {[30 30.05]}; {[40 40.05]}; {[50 50.05]}; {[60 60.05]}};           
% uses_tri={{a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}};

% Mawake330
% a=1:436;
% uses=   {{[1 1.05]}; {[2 2.05]}; {[4 4.05]}; {[6 6.05]}; {[8 8.05]}; {[10 10.05]}; {[12 12.05]}; {[14 14.05]}; {[16 16.05]}; {[18 18.05]}; {[20 20.05]}; {[30 30.05]}; {[40 40.05]}; {[50 50.05]}; {[60 60.05]}};           
% uses_tri={{a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}};

% Mawake378
% a=1:1360;
% uses=   {{[1 1.03]}; {[2 2.03]}; {[4 4.03]}; {[6 6.03]}; {[8 8.03]}; {[10 10.03]}; {[12 12.03]}; {[14 14.03]}; {[16 16.03]}; {[18 18.03]}; {[20 20.03]}; {[30 30.03]}; {[40 40.03]}; {[50 50.03]}; {[60 60.03]}};           
% uses_tri={{a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}};

% Mawake384
% a=1:1132;
% uses=   {{[0.01 5.01]}; {[0.02 5.02]}; {[0.04 5.04]}; {[0.06 5.06]}; {[0.08 5.08]}; {[0.10 5.1]}; {[0.12 5.12]}; {[0.14 5.14]}; {[0.16 5.16]}; {[0.18 5.18]}; {[0.20 5.20]}; {[0.30 5.30]}; {[0.40 5.40]}; {[0.50 5.50]}; {[0.60 5.60]}};           
% uses_tri={{a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}; {a a}};


% Mawake343
% uses=   {{1}; {2}; {3}};      % fileInds 7 to 66
% uses_tri={{[1:240]}; {[1:240]}; {[1:240]}};
% Mawake347
% uses=   {{1}; {2}; {3}};      % fileInds 20 to 109
% uses_tri={{[1:360]}; {[1:360]}; {[1:360]}};
% Mawake381
% uses=   {{1}; {2}; {3}};      % fileInds 6 to 66
% uses_tri={{[1:244]}; {[1:244]}; {[1:244]}};
% Mawake405
% uses=   {{1}};                % fileInds 2 to 89
% uses_tri={{[1:352]}};
% Mawake407
% uses=   {{1}};                % fileInds 2 to 63
% uses_tri={{[1:248]}};
% Mawake377
% uses=   {{1}; {2}; {3}};                % fileInds 10 to 73
% uses_tri={{[1:308]}; {[1:308]}; {[1:308]}};
% Mawake376
% uses=   {{1}; {2}; {3}};                % fileInds 10 to 73
% uses_tri={{[1:254]}; {[1:254]}; {[1:254]}};
% Mawake373
% uses=   {{1}; {2}; {3}};                % fileInds 25 to 79
% uses_tri={{[1:220]}; {[1:220]}; {[1:220]}};
% Mawake330
% uses=   {{1}; {2}};                % fileInds 3 to 60
% uses_tri={{[9:240]-8}; {[9:240]-8}};
% Mawake329
% uses=   {{1}; {2}};                % fileInds 201 to 223
% uses_tri={{[801:892]-800}; {[801:892]-800}};
% Mawake328
% uses=   {{1}; {2}};                % fileInds 189 to 217, 248 to 255
% uses_tri={{[753:868 989:1020]-752}; {[753:868 989:1020]-752}};
% Mawake327
% uses=   {{1}; {2}};                % fileInds 249 to 271
% uses_tri={{[993:1084]-992}; {[993:1084]-992}};
% Mawake344
% uses=   {{1}; {2}; {3}};                % fileInds 7 to 90, 97 to 102
% uses_tri={{[1:336 361:408]}; {[1:336 361:408]}; {[1:336 361:408]}};
% Mawake372
% uses=   {{1}; {2}; {3}};                % fileInds 1 to 32, 117 to 122
% uses_tri={{[1:128 465:488]}; {[1:128 465:488]}; {[1:128 465:488]}};
% Mawake375
% uses=   {{1}; {2}; {3}};                % fileInds 1 to 49
% uses_tri={{[1:196]}; {[1:196]}; {[1:196]}};
% Mawake342
% uses=   {{1}; {2}; {3}};                % fileInds 4 to 87
% uses_tri={{[1:336]}; {[1:336]}; {[1:336]}};
% Mawake362
% uses=   {{1}; {2}; {3}};                % fileInds 2 to 44, 82 to 98
% uses_tri={{[5:176 325:392]-4}; {[5:176 325:392]-4}; {[5:176 325:392]-4}};
% Mawake377
% uses=   {{1}; {2}; {3}};                % fileInds 1 to 73
% uses_tri={{[1:288]}; {[1:288]}; {[1:288]}};
% Mawake361
% uses=   {{1}; {2}; {3}};                             % fileInds 2 to 74
% uses_tri={{[1:292]}; {[1:292]}; {[1:292]}};
% Mawake354
% uses=   {{1}; {2}; {3}};                % fileInds 3 to 43, 99 to 106
% uses_tri={{[9:172 393:424]-8}; {[9:172 393:424]-8}; {[9:172 393:424]-8}};
% Mawake341
% uses=   {{1}; {2}; {3}};                              % fileInds 6 to 95
% uses_tri={{[21:380]-20}; {[21:380]-20}; {[21:380]-20}};
% Mawake353
% uses=   {{1}; {2}; {3}};                              % fileInds 1 to 93
% uses_tri={{[1:308]}; {[1:308]}; {[1:308]}};
% Mawake349
% uses=   {{1}; {2}; {3}};                              % fileInds 1 to 57
% uses_tri={{[1:228]}; {[1:228]}; {[1:228]}};
% Mawake355
% uses=   {{1}; {2}; {3}};                              % fileInds 6 to 92
% uses_tri={{[21:368]-20}; {[21:368]-20}; {[21:368]-20}};
% Mawake322
% uses=   {{1 1 1}; {1 1}};
% uses_tri={{[13:72]-12 [141:212]-12 [277:316]-12}; {[73:140]-12 [213:276]-12}};
% Mawake321
% uses=   {{1 1}; {1}};
% uses_tri={{[1:56] [157:200]}; {[57:156]}};
% Mawake320
% uses=   {{1 1}; {1 1}};
% uses_tri={{[1:80] [157:260]}; {[81:156] [261:340]}};
% Mawake319
% uses=   {{1 2 1}; {2 2}; {1}};
% uses_tri={{[1:80] [81:204] [273:332]}; {[1:80] [273:332]}; {[81:204]}};
% Mawake318
% uses=   {{1 1 1}; {2 2}; {2}};
% uses_tri={{[1:80] [81:196] [197:276]}; {[1:80] [197:276]}; {[81:196]}};
% Mawake317
% uses=   {{1 1}; {2 3}; {2}};
% uses_tri={{[1:188] [189:316 349:392]}; {[1:188] [189:316 349:392]}; {[189:316 349:392]}};
% Mawake316
% uses=   {{1 1}; {2}; {2}};
% uses_tri={{[1:136] [189:268]}; {[1:136]}; {[189:268]}};
% Mawake315
% uses=   {{1:2 1}; {3}; {2}};
% uses_tri={{[9:164]-8 [165:224]-8}; {[9:164]-8}; {[165:224]-8}};
% Mawake314
% uses=   {{1:2 1 1}; {3}; {2}; {2}};
% uses_tri={{[1:160] [197:268] [269:308]}; {[1:160]}; {[197:268]}; {[269:308]}};
% Mawake313
% uses=   {{1:2 1}; {3}; {2}};
% uses_tri={{[1:128] [169:244]}; {[1:128]}; {[169:244]}};
% Mawake312
% uses=   {{1:2 1}; {3}; {2}};
% uses_tri={{[33:132]-32 [169:232]-32}; {[33:132]-32}; {[169:232]-32}};
% Mawake311
% uses=   {{1:2 1 1}; {3}; {2}; {2}};
% uses_tri={{[13:148]-12 [204:254]-12 [255:299]-12}; {[13:148]-12}; {[204:254]-12}; {[255:299]-12}};
% Mawake310
% uses=   {{1:2 1:2 1}; {3 3}; {2}};
% uses_tri={{[9:144]-8 [201:228]-8 [229:279]-8}; {[9:144]-8 [201:228]-8}; {[229:279]-8}};
% Mawake206
% uses=   {{1}; {2}; {3}; {4}; {5}; {6}; {7}; {8}; {9}; {10}; {11}; {12}};
% uses_tri={{[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}; {[1:252]}};


%uses=unique(dLGNpsth.unitStimcond{1});
% usealls=1:10000;
usealls=[1 1.05 2 2.05 4 4.05 6 6.05 8 8.05 10 10.05 12 12.05 14 14.05 16 16.05 18 18.05 20 20.05 30 30.05 40 40.05 50 50.05 60 60.05];

% usealls=[1.0000 ...
%     1.0100 ...
%     1.0300 ...
%     2.0000 ...
%     2.0100 ...
%     2.0300 ...
%     4.0000 ...
%     4.0100 ...
%     4.0300 ...
%     6.0000 ...
%     6.0100 ...
%     6.0300 ...
%     8.0000 ...
%     8.0100 ...
%     8.0300 ...
%    10.0000 ...
%    10.0100 ...
%    10.0300 ...
%    12.0000 ...
%    12.0100 ...
%    12.0300 ...
%    14.0000 ...
%    14.0100 ...
%    14.0300 ...
%    16.0000 ...
%    16.0100 ...
%    16.0300 ...
%    18.0000 ...
%    18.0100 ...
%    18.0300 ...
%    20.0000 ...
%    20.0100 ...
%    20.0300 ...
%    30.0000 ...
%    30.0100 ...
%    30.0300 ...
%    40.0000 ...
%    40.0100 ...
%    40.0300 ...
%    50.0000 ...
%    50.0100 ...
%    50.0300 ...
%    60.0000 ...
%    60.0100 ...
%    60.0300];

% usealls=[0.0100 ...
%     0.0200 ...
%     0.0400 ...
%     0.0600 ...
%     0.0800 ...
%     0.1000 ...
%     0.1200 ...
%     0.1400 ...
%     0.1600 ...
%     0.1800 ...
%     0.2000 ...
%     0.3000 ...
%     0.4000 ...
%     0.5000 ...
%     0.6000 ...
%     5.0100 ...
%     5.0200 ...
%     5.0400 ...
%     5.0600 ...
%     5.0800 ...
%     5.1000 ...
%     5.1200 ...
%     5.1400 ...
%     5.1600 ...
%     5.1800 ...
%     5.2000 ...
%     5.3000 ...
%     5.4000 ...
%     5.5000 ...
%     5.6000];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if isempty(LFPdata)
%     if isempty(LFP_Fs)
%         LFP_Fs=2500;
%     end
%     
%     disp('1');
%     [~,~,LFPbySweep,ledConds,stimConds]=makeCSD2(expt,useFileInd,[],[],[]);
%     save([outputDir '\' 'LFPbySweep.mat'],'LFPbySweep'); % V1 then hippo
%     save([outputDir '\' 'LFPbySweep_ledConds.mat'],'ledConds');
%     save([outputDir '\' 'LFPbySweep_stimConds.mat'],'stimConds');
%     
%     L=LFPbySweep{1};
%     L=L(~isnan(ledConds) & ~isnan(stimConds),:);
%     LFPbySweep{1}=L;
%     L=LFPbySweep{2};
%     L=L(~isnan(ledConds) & ~isnan(stimConds),:);
%     LFPbySweep{2}=L;
% else
%     LFPbySweep{1}=LFPdata.LFPbySweep{1};
%     if length(LFPdata.LFPbySweep)>1
%         LFPbySweep{2}=LFPdata.LFPbySweep{2};
%     end
%     ledConds=LFPdata.ledConds;
%     stimConds=LFPdata.stimConds;
%     
%     L=LFPbySweep{1};
%     L=L(~isnan(ledConds) & ~isnan(stimConds),:);
%     LFPbySweep{1}=L;
%     if length(LFPdata.LFPbySweep)>1
%         L=LFPbySweep{2};
%         L=L(~isnan(ledConds) & ~isnan(stimConds),:);
%         LFPbySweep{2}=L;
%     end
% end
% 
% if size(LFPbySweep{1},1)~=length(unique(dLGNpsth.unitTrials{1}))
%     disp('Mismatched number of trials in LFPbySweep and dLGNpsth');
%     return
% end
% 
% disp('2');
% if length(LFPbySweep)>1
%     hippoLFP{1}=LFPbySweep{2};
% else
%     hippoLFP{1}=LFPbySweep{1};
% end
% thetaDiff=getHippoTheta(hippoLFP,ledConds,stimConds,usel_noLED,usealls,LFP_Fs);
% figure(); 
% hist(nanmean(thetaDiff,2),30);
% title('Histogram of theta diff');
% 
% cxSynch=[];
% if length(LFPbySweep)>1
%     cxLFP{1}=LFPbySweep{1};
%     cxSynch=getCortexSynch(cxLFP,ledConds,stimConds,usel_noLED,usealls,LFP_Fs);
% end
% 
% if length(LFPbySweep)>1
%     figure(); 
%     scatter(nanmean(thetaDiff,2),nanmean(cxSynch,2)); 
% end
% save([outputDir '\' 'thetaDiff_noLED.mat'],'thetaDiff');
% save([outputDir '\' 'cxSynch_noLED.mat'],'cxSynch');
% save([outputDir '\' 'dLGNpsth.mat'],'dLGNpsth');
% 
% thetaDiff=getHippoTheta(hippoLFP,ledConds,stimConds,[usel_noLED usel_LED],usealls,LFP_Fs);
% save([outputDir '\' 'thetaDiff.mat'],'thetaDiff');
% 
% thetaThresh=input('Theta difference threshold to use?');
% save([outputDir '\' 'thetaThresh.mat'],'thetaThresh');
% 
% noThetaTrials=nanmean(thetaDiff,2)<thetaThresh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([outputDir '\' 'noThetaTrials.mat'],'noThetaTrials');
save([outputDir '\' 'dLGNpsth.mat'],'dLGNpsth');

disp('3');
tri=dLGNpsth.unitTrials{1};
if length(noThetaTrials)~=length(tri)
    disp('Length of noThetaTrials does not match number of trials in dLGNpsth');
    return
end
filt_psth=filtPSTH(dLGNpsth,tri(noThetaTrials'));
filt_psth_theta=filtPSTH(dLGNpsth,tri(noThetaTrials'~=1));

% No LED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usel=usel_noLED;

% uses=   {{1 1 1}; {1 1}};
% uses_tri={{[13:72]-12 [141:212]-12 [277:316]-12}; {[73:140]-12 [213:276]-12}};

for i=1:length(uses)
    currs=uses{i};
    currtri=uses_tri{i};
    allcurr_psth=[];
    allcurr_psth_theta=[];
    for j=1:length(currs)
        subcurrs=currs{j};
        subcurrtri=currtri{j};
        
        sforsub=filt_psth.unitStimcond{1};
        triforsub=filt_psth.unitTrials{1};
        curr_psth=filtPSTH(filt_psth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        sforsub=filt_psth_theta.unitStimcond{1};
        triforsub=filt_psth_theta.unitTrials{1};
        curr_psth_theta=filtPSTH(filt_psth_theta,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        allcurr_psth=concatPSTHs(allcurr_psth,curr_psth);
        allcurr_psth_theta=concatPSTHs(allcurr_psth_theta,curr_psth_theta);
    end    
    [noTheta_trialAv_temp(i).allS,noTheta_trialAv_temp(i).HFa,noTheta_trialAv_temp(i).LFa,noTheta_trialAv_temp(i).F1amp,noTheta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),1,usealls,usel,freqs(i));
    [theta_trialAv_temp(i).allS,theta_trialAv_temp(i).HFa,theta_trialAv_temp(i).LFa,theta_trialAv_temp(i).F1amp,theta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),1,usealls,usel,freqs(i));
end
save([outputDir '\' 'noTheta_trialAv_temp_noLED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_noLED.mat'],'theta_trialAv_temp');
noTheta_trialAv=avResponses(noTheta_trialAv_temp);
theta_trialAv=avResponses(theta_trialAv_temp);
save([outputDir '\' 'noTheta_trialAv_noLED.mat'],'noTheta_trialAv');
save([outputDir '\' 'theta_trialAv_noLED.mat'],'theta_trialAv');

[noTheta.allS,noTheta.HFa,noTheta.LFa,noTheta.F1amp,noTheta.allpower]=getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel,1);
save([outputDir '\' 'noTheta_noLED.mat'],'noTheta');
[theta.allS,theta.HFa,theta.LFa,theta.F1amp,theta.allpower]=getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel,1);
save([outputDir '\' 'theta_noLED.mat'],'theta');

disp('4');
% With LED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usel=usel_LED;
for i=1:length(uses)
    currs=uses{i};
    currtri=uses_tri{i};
    allcurr_psth=[];
    allcurr_psth_theta=[];
    for j=1:length(currs)
        subcurrs=currs{j};
        subcurrtri=currtri{j};
        
        sforsub=filt_psth.unitStimcond{1};
        triforsub=filt_psth.unitTrials{1};
        curr_psth=filtPSTH(filt_psth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        sforsub=filt_psth_theta.unitStimcond{1};
        triforsub=filt_psth_theta.unitTrials{1};
        curr_psth_theta=filtPSTH(filt_psth_theta,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        allcurr_psth=concatPSTHs(allcurr_psth,curr_psth);
        allcurr_psth_theta=concatPSTHs(allcurr_psth_theta,curr_psth_theta);
    end    
    [noTheta_trialAv_temp(i).allS,noTheta_trialAv_temp(i).HFa,noTheta_trialAv_temp(i).LFa,noTheta_trialAv_temp(i).F1amp,noTheta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),1,usealls,usel,freqs(i));
    [theta_trialAv_temp(i).allS,theta_trialAv_temp(i).HFa,theta_trialAv_temp(i).LFa,theta_trialAv_temp(i).F1amp,theta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth_theta,ones(length(allcurr_psth_theta.psths),1),zeros(length(allcurr_psth_theta.psths),1),1,usealls,usel,freqs(i));
end
save([outputDir '\' 'noTheta_trialAv_temp_LED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_LED.mat'],'theta_trialAv_temp');
noTheta_trialAv=avResponses(noTheta_trialAv_temp);
theta_trialAv=avResponses(theta_trialAv_temp);
save([outputDir '\' 'noTheta_trialAv_LED.mat'],'noTheta_trialAv');
save([outputDir '\' 'theta_trialAv_LED.mat'],'theta_trialAv');

[noTheta.allS,noTheta.HFa,noTheta.LFa,noTheta.F1amp,noTheta.allpower]=getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel,1);
save([outputDir '\' 'noTheta_LED.mat'],'noTheta');
[theta.allS,theta.HFa,theta.LFa,theta.F1amp,theta.allpower]=getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel,1);
save([outputDir '\' 'theta_LED.mat'],'theta');

end

function outresponse=avResponses(response)

HFa=zeros(size(response(1).HFa));
LFa=zeros(size(response(1).LFa));
F1amp=zeros(size(response(1).F1amp));
allpower=zeros(size(response(1).allpower));
for i=1:length(response)
    if ~isnan(response(i).HFa)
        HFa=HFa+response(i).HFa;
    end
    if ~isnan(response(i).LFa)
        LFa=LFa+response(i).LFa;
    end
    if ~isnan(response(i).F1amp)
        F1amp=F1amp+response(i).F1amp;
    end
    if ~isnan(response(i).allpower)
        allpower=allpower+response(i).allpower;
    end
end
HFa=HFa./length(response);
LFa=LFa./length(response);
F1amp=F1amp./length(response);
allpower=allpower./length(response);
outresponse.HFa=HFa;
outresponse.LFa=LFa;
outresponse.F1amp=F1amp;
outresponse.allpower=allpower;

end


function newPSTH=concatPSTHs(struct1,struct2)

if isempty(struct1)
    newPSTH=struct2;
    return
elseif isempty(struct2)
    newPSTH=struct1;
    return
else
    if length(struct1.t)~=length(struct2.t)
        disp('t in psths does not match');
        return
    end
    newPSTH.t=struct1.t;
    newPSTH.psths=cell(length(struct1.psths),1);
    for i=1:length(struct1.psths)
        newPSTH.psths{i}=[struct1.psths{i}; struct2.psths{i}];
        newPSTH.unitTrials{i}=[struct1.unitTrials{i} struct2.unitTrials{i}];
        newPSTH.unitStimcond{i}=[struct1.unitStimcond{i} struct2.unitStimcond{i}];
        newPSTH.unitLED{i}=[struct1.unitLED{i} struct2.unitLED{i}];
    end
end

end


