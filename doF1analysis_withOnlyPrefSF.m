function doF1analysis_withOnlyPrefSF(expt,useFileInd,outputDir,LFP_Fs,dLGNpsth,usel_noLED,usel_LED,LFPdata,uses,uses_tri,noThetaTrials)

windowbeforestim=[0 3];
windowduringstim=[4.5 6.5];

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
% uses_tri={{[5:176 325:392]}; {[5:176 325:392]}; {[5:176 325:392]}};
% Mawake377
% uses=   {{1}; {2}; {3}};
% uses_tri={{[1:288]}; {[1:288]}; {[1:288]}};
% Mawake361
% uses=   {{1}; {2}; {3}};                             % fileInds 2 to 74
% uses_tri={{[1:292]}; {[1:292]}; {[1:292]}};
% Mawake354
% uses=   {{1}; {2}; {3}};                % fileInds 3 to 43, 99 to 106
% uses_tri={{[9:172 393:424]}; {[9:172 393:424]}; {[9:172 393:424]}};
% Mawake341
% uses=   {{1}; {2}; {3}};                              % fileInds 6 to 95
% uses_tri={{[21:380]}; {[21:380]}; {[21:380]}};
% Mawake353
% uses=   {{1}; {2}; {3}};                              % fileInds 1 to 93
% uses_tri={{[1:308]}; {[1:308]}; {[1:308]}};
% Mawake349
% uses=   {{1}; {2}; {3}};                              % fileInds 1 to 57
% uses_tri={{[1:228]}; {[1:228]}; {[1:228]}};
% Mawake355
% uses=   {{1}; {2}; {3}};                              % fileInds 6 to 92
% uses_tri={{[21:368]}; {[21:368]}; {[21:368]}};
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
% uses_tri={{[9:164] [165:224]}; {[9:164]}; {[165:224]}};
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
usealls=1:10000;

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

% ADDED TO PICK PREF RESPONSE
visresponse=nan(length(uses),length(filt_psth.psths));
% END ADDED TO PICK PREF RESPONSE
for i=1:length(uses)
    currs=uses{i};
    currtri=uses_tri{i};
    allcurr_psth=[];
    allcurr_psth_theta=[];
    all_forvisresponse=[];
    for j=1:length(currs)
        subcurrs=currs{j};
        subcurrtri=currtri{j};
        
        sforsub=filt_psth.unitStimcond{1};
        triforsub=filt_psth.unitTrials{1};
        curr_psth=filtPSTH(filt_psth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        sforsub=filt_psth_theta.unitStimcond{1};
        triforsub=filt_psth_theta.unitTrials{1};
        curr_psth_theta=filtPSTH(filt_psth_theta,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        
        % ADDED TO PICK PREF RESPONSE
        sforsub=dLGNpsth.unitStimcond{1};
        triforsub=dLGNpsth.unitTrials{1};
        curr_psth_for_visresponse=filtPSTH(dLGNpsth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        for psthcount=1:length(curr_psth_for_visresponse.psths)
            currpsthcount=curr_psth_for_visresponse.psths{psthcount};
            visresponse(i,psthcount)=nanmean(nanmean(currpsthcount(:,dLGNpsth.t>=windowduringstim(1) & dLGNpsth.t<=windowduringstim(2)),2),1)-nanmean(nanmean(currpsthcount(:,dLGNpsth.t>=windowbeforestim(1) & dLGNpsth.t<=windowbeforestim(2)),2),1);
        end
        % END ADDED TO PICK PREF RESPONSE
        
        allcurr_psth=concatPSTHs(allcurr_psth,curr_psth);
        allcurr_psth_theta=concatPSTHs(allcurr_psth_theta,curr_psth_theta);
        all_forvisresponse=concatPSTHs(all_forvisresponse,curr_psth_for_visresponse);
    end    
    
%     % ADDED TO PICK PREF RESPONSE
%     [~,ismaxstim]=max(visresponse,[],1);
%     allcurr_psth.t=filt_psth.t;
%     allcurr_psth_theta.t=filt_psth_theta.t;
%     for j=1:length(currs)
%         subcurrs=currs{j};
%         subcurrtri=currtri{j};
%         
%         sforsub=filt_psth.unitStimcond{1};
%         triforsub=filt_psth.unitTrials{1};
%         curr_psth=filtPSTH(filt_psth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
%         
%         sforsub=filt_psth_theta.unitStimcond{1};
%         triforsub=filt_psth_theta.unitTrials{1};
%         curr_psth_theta=filtPSTH(filt_psth_theta,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
%         
%         for psthcount=1:length(filt_psth.psths)
%             if ismaxstim(psthcount)==j
%                 allcurr_psth.psths{psthcount}=curr_psth.psths{psthcount};
%                 allcurr_psth.unitTrials{psthcount}=curr_psth.unitTrials{psthcount};
%                 allcurr_psth.unitStimcond{psthcount}=curr_psth.unitStimcond{psthcount};
%                 allcurr_psth.unitLED{psthcount}=curr_psth.unitLED{psthcount};
%                 
%                 allcurr_psth_theta.psths{psthcount}=curr_psth_theta.psths{psthcount};
%                 allcurr_psth_theta.unitTrials{psthcount}=curr_psth_theta.unitTrials{psthcount};
%                 allcurr_psth_theta.unitStimcond{psthcount}=curr_psth_theta.unitStimcond{psthcount};
%                 allcurr_psth_theta.unitLED{psthcount}=curr_psth_theta.unitLED{psthcount};
%             end
%         end
%     end
%     % END ADDED TO PICK PREF RESPONSE         
        
        
        
    [noTheta_trialAv_temp(i).allS,noTheta_trialAv_temp(i).HFa,noTheta_trialAv_temp(i).LFa,noTheta_trialAv_temp(i).F1amp,noTheta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),1,usealls,usel);
    [theta_trialAv_temp(i).allS,theta_trialAv_temp(i).HFa,theta_trialAv_temp(i).LFa,theta_trialAv_temp(i).F1amp,theta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),1,usealls,usel);
    [~,~,~,visresponseF1(i).F1amp,~]=getHFandLFalphaResponses(all_forvisresponse,ones(length(all_forvisresponse.psths),1),zeros(length(all_forvisresponse.psths),1),1,usealls,usel);
end
save([outputDir '\' 'noTheta_trialAv_temp_noLED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_noLED.mat'],'theta_trialAv_temp');
save([outputDir '\' 'visresponse_noLED.mat'],'visresponse');

noTheta_trialAv=pickPrefResponse(noTheta_trialAv_temp,visresponse);
theta_trialAv=pickPrefResponse(theta_trialAv_temp,visresponse);
% noTheta_trialAv=avResponses(noTheta_trialAv_temp);
% theta_trialAv=avResponses(theta_trialAv_temp);
save([outputDir '\' 'noTheta_trialAv_noLED.mat'],'noTheta_trialAv');
save([outputDir '\' 'theta_trialAv_noLED.mat'],'theta_trialAv');

[noTheta.allS,noTheta.HFa,noTheta.LFa,noTheta.F1amp,noTheta.allpower]=getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel);
save([outputDir '\' 'noTheta_noLED.mat'],'noTheta');
[theta.allS,theta.HFa,theta.LFa,theta.F1amp,theta.allpower]=getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel);
save([outputDir '\' 'theta_noLED.mat'],'theta');

disp('4');
% With LED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usel=usel_LED;
% ADDED TO PICK PREF RESPONSE
visresponse=nan(length(uses),length(filt_psth.psths));
% END ADDED TO PICK PREF RESPONSE
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
        
        % ADDED TO PICK PREF RESPONSE
        sforsub=dLGNpsth.unitStimcond{1};
        triforsub=dLGNpsth.unitTrials{1};
        curr_psth_for_visresponse=filtPSTH(dLGNpsth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
        for psthcount=1:length(curr_psth_for_visresponse.psths)
            currpsthcount=curr_psth_for_visresponse.psths{psthcount};
            visresponse(i,psthcount)=nanmean(nanmean(currpsthcount(:,dLGNpsth.t>=windowduringstim(1) & dLGNpsth.t<=windowduringstim(2)),2),1)-nanmean(nanmean(currpsthcount(:,dLGNpsth.t>=windowbeforestim(1) & dLGNpsth.t<=windowbeforestim(2)),2),1);
        end
        % END ADDED TO PICK PREF RESPONSE
        
        allcurr_psth=concatPSTHs(allcurr_psth,curr_psth);
        allcurr_psth_theta=concatPSTHs(allcurr_psth_theta,curr_psth_theta);
    end
    
    
    
%     % ADDED TO PICK PREF RESPONSE
%     [~,ismaxstim]=max(visresponse,[],1);
%     allcurr_psth.t=filt_psth.t;
%     allcurr_psth_theta.t=filt_psth_theta.t;
%     for j=1:length(currs)
%         subcurrs=currs{j};
%         subcurrtri=currtri{j};
%         
%         sforsub=filt_psth.unitStimcond{1};
%         triforsub=filt_psth.unitTrials{1};
%         curr_psth=filtPSTH(filt_psth,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
%         
%         sforsub=filt_psth_theta.unitStimcond{1};
%         triforsub=filt_psth_theta.unitTrials{1};
%         curr_psth_theta=filtPSTH(filt_psth_theta,triforsub(ismember(triforsub,subcurrtri) & ismember(sforsub,subcurrs)));
%         
%         for psthcount=1:length(filt_psth.psths)
%             if ismaxstim(psthcount)==j
%                 allcurr_psth.psths{psthcount}=curr_psth.psths{psthcount};
%                 allcurr_psth.unitTrials{psthcount}=curr_psth.unitTrials{psthcount};
%                 allcurr_psth.unitStimcond{psthcount}=curr_psth.unitStimcond{psthcount};
%                 allcurr_psth.unitLED{psthcount}=curr_psth.unitLED{psthcount};
%                 
%                 allcurr_psth_theta.psths{psthcount}=curr_psth_theta.psths{psthcount};
%                 allcurr_psth_theta.unitTrials{psthcount}=curr_psth_theta.unitTrials{psthcount};
%                 allcurr_psth_theta.unitStimcond{psthcount}=curr_psth_theta.unitStimcond{psthcount};
%                 allcurr_psth_theta.unitLED{psthcount}=curr_psth_theta.unitLED{psthcount};
%             end
%         end
%     end
%     % END ADDED TO PICK PREF RESPONSE   
    
    
    
    [noTheta_trialAv_temp(i).allS,noTheta_trialAv_temp(i).HFa,noTheta_trialAv_temp(i).LFa,noTheta_trialAv_temp(i).F1amp,noTheta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),1,usealls,usel);
    [theta_trialAv_temp(i).allS,theta_trialAv_temp(i).HFa,theta_trialAv_temp(i).LFa,theta_trialAv_temp(i).F1amp,theta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth_theta,ones(length(allcurr_psth_theta.psths),1),zeros(length(allcurr_psth_theta.psths),1),1,usealls,usel);
end
save([outputDir '\' 'noTheta_trialAv_temp_LED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_LED.mat'],'theta_trialAv_temp');
save([outputDir '\' 'visresponse_LED.mat'],'visresponse');

noTheta_trialAv=pickPrefResponse(noTheta_trialAv_temp,visresponse);
theta_trialAv=pickPrefResponse(theta_trialAv_temp,visresponse);
% noTheta_trialAv=avResponses(noTheta_trialAv_temp);
% theta_trialAv=avResponses(theta_trialAv_temp);
save([outputDir '\' 'noTheta_trialAv_LED.mat'],'noTheta_trialAv');
save([outputDir '\' 'theta_trialAv_LED.mat'],'theta_trialAv');

[noTheta.allS,noTheta.HFa,noTheta.LFa,noTheta.F1amp,noTheta.allpower]=getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel);
save([outputDir '\' 'noTheta_LED.mat'],'noTheta');
[theta.allS,theta.HFa,theta.LFa,theta.F1amp,theta.allpower]=getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel);
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

function outresponse=pickPrefResponse(response,visresponse)

[~,ismaxstim]=nanmax(visresponse,[],1);
% [~,ismaxstim]=nanmin(visresponse,[],1);
outresponse.HFa=nan(size(response(1).HFa));
outresponse.LFa=nan(size(response(1).LFa));
outresponse.F1amp=nan(size(response(1).F1amp));
outresponse.allpower=nan(size(response(1).allpower));
for i=1:length(ismaxstim)
    currmaxstim=ismaxstim(i);
    outresponse.HFa(i,:)=response(currmaxstim).HFa(i,:);
    outresponse.LFa(i,:)=response(currmaxstim).LFa(i,:);
    outresponse.F1amp(i,:)=response(currmaxstim).F1amp(i,:);
    outresponse.allpower(i,:)=response(currmaxstim).allpower(i,:);
end
   

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


