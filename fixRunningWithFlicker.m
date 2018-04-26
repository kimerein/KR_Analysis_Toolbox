function [velocity,spikes]=fixRunningWithFlicker(runningData,ledData,spikes,velocity)

trialDuration=3.5;
downSampFactor=250;
% scaleFactor=0.9392;
discontinuityThresh=0.2;
% runningThresh=0.02;
runningThresh=0.05;

if isempty(velocity)
    runbase=mean(runningData(1,1:2000));
    ledbase=mean(ledData(1,1:2000));
    runningData=runningData-runbase;
    ledData=ledData-ledbase;
    maRun=max(runningData(1,:));
    maLed=max(ledData(1,:));
    runningData=runningData-ledData.*1.0.*(maRun/maLed);
    
    % Patch for m116
    runningData(:,3631:end)=runningData(:,3631:end)-repmat(runningData(:,3631)-runningData(:,3617),1,size(runningData(:,3631:end),2));
    
    % Determine animal's velocity
    r=runningData;
    
    r(:,1:10)=repmat(r(:,11),1,10);
    % for i=1:size(r,1)
    %     r(i,1:10)=r(i,11)*ones(1,10);
    % end
    
    % Remove discontinuities in position data
    rdiff=diff(r');
    rdiff=rdiff';
    for i=1:size(rdiff,1)
        inds=find(abs(rdiff(i,:))>=discontinuityThresh);
        useinds=inds(inds>10);
        for j=useinds
            if j+4>size(r,2)
                r(i,end-4:end)=r(i,end-4)*ones(1,5);
            else
                r(i,:)=[r(i,1:j) r(i,j) r(i,j) r(i,j) r(i,j+4:size(r,2))-r(i,j+4)+r(i,j)];
            end
        end
    end
    
    % Down sample running signal because animal can't run above a certain speed
    r=downSampMatrix(r,downSampFactor);
    
    der=r(:,1:end-floor(1000/downSampFactor))-r(:,floor(1000/downSampFactor)+1:end); % 1000 fewer points than r per row
    velocity=[der repmat(der(:,size(der,2)),1,floor(1000/downSampFactor))];
    % for i=1:size(der,1)
    %     velocity(i,:)=[der(i,:) der(i,end)*ones(1,1000)];
    % end
    velocity=abs(velocity);
    % for i=1:size(velocity,1)
    %     velocity(i,:)=smooth(velocity(i,:)); % default window is 5
    % end
    
    
    return
else
    % Make field in spikes that identifies running
    t=unique(spikes.sweeps.trials);
    if length(t)~=size(velocity,1)
        disp('trials and velocity do not match');
        return
    end
    
    t=sort(t);
    spikes.running=nan(size(spikes.spiketimes));
    v_x=linspace(0,trialDuration,size(velocity,2));
    for i=1:length(t)
        st=spikes.spiketimes(spikes.trials==t(i));
        isrun=nan(1,length(st));
        for j=1:length(st)
            f=find(v_x>=st(j),1,'first');
            isrun(j)=velocity(i,f)>=runningThresh;
        end
        spikes.running(spikes.trials==t(i))=isrun;
    end
end
