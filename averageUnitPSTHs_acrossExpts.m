function [psth,blackPSTH,redPSTH]=averageUnitPSTHs_acrossExpts(dd,psthName,infileName,tt,thetaStateName,infile_thetaStateName,take)

useThetaState=0;

psth=[];
for i=1:length(dd)
    d=dd{i};
    if ~exist([d '\' psthName '.mat'],'file')
        disp([d '\' psthName ' does not exist ... Skipping']);
        continue
    end
    a=load([d '\' psthName]);
    currpsth=a.(infileName);
    if exist([d '\' take '.mat'],'file')
        a=load([d '\' take '.mat']);
        if ~isfield(a,'take')
            error('variable must be called "take" in file');
        else
            takeTri=a.take;
        end
    else
        takeTri=[];
    end
    if strcmp(infileName,'fractionpsth')
        currpsth.psths=currpsth.psths';
    end
    if useThetaState==1
        a=load([tt{i} '\' thetaStateName]);
        currtheta=a.(infile_thetaStateName);
        if length(currtheta)>size(currpsth.psths{1},1) % this is a hack
            currtheta=currtheta(1:size(currpsth.psths{1},1));
        end
        for j=1:length(currpsth.psths)
            currpsth.noThetaTrials{j}=currtheta;
        end
    end
    
    % take only certain trials?
    if ~isempty(takeTri)
        for j=1:length(currpsth.psths)
            temp=currpsth.psths{j};
            temp=temp(takeTri,:);
            currpsth.psths{j}=temp;
            
            temp=currpsth.unitTrials{j};
            temp=temp(takeTri);
            currpsth.unitTrials{j}=temp;
            
            temp=currpsth.unitStimcond{j};
            temp=temp(takeTri);
            currpsth.unitStimcond{j}=temp;
            
            temp=currpsth.unitLED{j};
            temp=temp(takeTri);
            currpsth.unitLED{j}=temp;
            
            temp=currpsth.noThetaTrials{j};
            temp=temp(takeTri);
            currpsth.noThetaTrials{j}=temp;
        end
    end
    
    psth=concatPSTHs_acrossExpts(psth,currpsth);
end

% take=[9    10    13    24    31    34    45    48    52    59    96    97   121 ...
%    125   133   134   138   145   151   153   161   162   168   169   170   172  174];

% psth.psths=psth.psths(take);
% psth.unitTrials=psth.unitTrials(take);
% psth.unitStimcond=psth.unitStimcond(take);
% psth.unitLED=psth.unitLED(take);

if useThetaState==1
    noTheta_psth.t=psth.t;
    for i=1:length(psth.psths)
        thetaState=psth.noThetaTrials{i};
        p=psth.psths{i};
        noTheta_psth.psths{i}=p(thetaState==1,:);
        p=psth.unitTrials{i};
        noTheta_psth.unitTrials{i}=p(thetaState==1);
        p=psth.unitStimcond{i};
        noTheta_psth.unitStimcond{i}=p(thetaState==1);
        p=psth.unitLED{i};
        noTheta_psth.unitLED{i}=p(thetaState==1);
    end
    
    theta_psth.t=psth.t;
    for i=1:length(psth.psths)
        thetaState=psth.noThetaTrials{i};
        p=psth.psths{i};
        theta_psth.psths{i}=p(thetaState==0,:);
        p=psth.unitTrials{i};
        theta_psth.unitTrials{i}=p(thetaState==0);
        p=psth.unitStimcond{i};
        theta_psth.unitStimcond{i}=p(thetaState==0);
        p=psth.unitLED{i};
        theta_psth.unitLED{i}=p(thetaState==0);
    end
    [blackPSTH,redPSTH]=averageUnitPSTHs(noTheta_psth);
%     [blackPSTH,redPSTH]=averageUnitPSTHs(theta_psth);
else
    [blackPSTH,redPSTH]=averageUnitPSTHs(psth);
end

end

function newPSTH=concatPSTHs_acrossExpts(struct1,struct2)

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
    newPSTH.psths=[struct1.psths; struct2.psths];
    if size(struct2.unitTrials,1)>1
        struct2.unitTrials=struct2.unitTrials';
        struct2.unitStimcond=struct2.unitStimcond';
        struct2.unitLED=struct2.unitLED';
    end
    newPSTH.unitTrials=[struct1.unitTrials struct2.unitTrials];
    newPSTH.unitStimcond=[struct1.unitStimcond struct2.unitStimcond];
    newPSTH.unitLED=[struct1.unitLED struct2.unitLED];
    if isfield(struct1, 'noThetaTrials')
        newPSTH.noThetaTrials=[struct1.noThetaTrials struct2.noThetaTrials];
    end
end

end