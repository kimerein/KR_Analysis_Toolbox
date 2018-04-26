function [xpoints,ypoints,grp_result]=calcUnitPersistence(spikes,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED)

% grp_led{1}=[4 8];
% grp_window{1}=[1.34 1.5]; % d
% grp_stimcond{1}=1:24;
% grp_led{2}=[4 8];
% grp_window{2}=[1.18 1.34]; % b
% grp_stimcond{2}=1:24;
% grp_led{3}=[3 5];
% grp_window{3}=[1.34 1.5]; % c
% grp_stimcond{3}=1:24;
% grp_led{4}=[3 5];
% grp_window{4}=[1.18 1.34]; % a
% grp_stimcond{4}=1:24;
% grp_led{1}=[5];
% grp_window{1}=[1.3 1.46]; % d
% grp_stimcond{1}=[101];
% grp_led{2}=[5];
% grp_window{2}=[1.14 1.3]; % b
% grp_stimcond{2}=[101];
% grp_led{3}=[5];
% grp_window{3}=[1.3 1.46]; % c
% grp_stimcond{3}=[100];
% grp_led{4}=[5];
% grp_window{4}=[1.14 1.3]; % a
% grp_stimcond{4}=[100];

grp_led{1}=[4 8];
grp_window{1}=[1.34 1.5]; % d
grp_stimcond{1}=1:8;
grp_led{2}=[4 8];
grp_window{2}=[1.18 1.34]; % b
grp_stimcond{2}=1:8;
grp_led{3}=[3 5];
grp_window{3}=[1.34 1.5]; % c
grp_stimcond{3}=1:8;
grp_led{4}=[3 5];
grp_window{4}=[1.18 1.34]; % a
grp_stimcond{4}=1:8;
grp_led{5}=[4 8];
grp_window{5}=[0.6 1]; % base_b
grp_stimcond{5}=1:8;
grp_led{6}=[3 5];
grp_window{6}=[0.6 1]; % base_a
grp_stimcond{6}=1:8;
% grp_led{1}=[5];
% grp_window{1}=[1.3 1.46]; % d
% grp_stimcond{1}=[101];
% grp_led{2}=[5];
% grp_window{2}=[1.14 1.3]; % b
% grp_stimcond{2}=[101];
% grp_led{3}=[5];
% grp_window{3}=[1.3 1.46]; % c
% grp_stimcond{3}=[100];
% grp_led{4}=[5];
% grp_window{4}=[1.14 1.3]; % a
% grp_stimcond{4}=[100];
% grp_led{5}=[5];
% grp_window{5}=[0.6 1]; % base_b
% grp_stimcond{5}=[101];
% grp_led{6}=[5];
% grp_window{6}=[0.6 1]; % base_a
% grp_stimcond{6}=[100];

% grp_led{1}=[4 8];
% grp_window{1}=[1.08 1.111]; % d
% grp_stimcond{1}=1:8;
% grp_led{2}=[4 8];
% grp_window{2}=[1.111 1.142]; % b
% grp_stimcond{2}=1:8;
% grp_led{3}=[3 5];
% grp_window{3}=[1.08 1.111]; % c
% grp_stimcond{3}=1:8;
% grp_led{4}=[3 5];
% grp_window{4}=[1.111 1.142]; % a
% grp_stimcond{4}=1:8;
% grp_led{5}=[4 8];
% grp_window{5}=[0.6 1]; % base_b
% grp_stimcond{5}=1:8;
% grp_led{6}=[3 5];
% grp_window{6}=[0.6 1]; % base_a
% grp_stimcond{6}=1:8;
% grp_led{1}=[4 8];
% grp_window{1}=[1.06 1.2]; % d
% grp_stimcond{1}=1:8;
% grp_led{2}=[4 8];
% grp_window{2}=[1.2 1.34]; % b
% grp_stimcond{2}=1:8;
% grp_led{3}=[3 5];
% grp_window{3}=[1.06 1.2]; % c
% grp_stimcond{3}=1:8;
% grp_led{4}=[3 5];
% grp_window{4}=[1.2 1.34]; % a
% grp_stimcond{4}=1:8;
% grp_led{5}=[4 8];
% grp_window{5}=[0.6 1]; % base_b
% grp_stimcond{5}=1:8;
% grp_led{6}=[3 5];
% grp_window{6}=[0.6 1]; % base_a
% grp_stimcond{6}=1:8;

% persistence_y=@(a,b,c,d) d./b;
% persistence_x=@(a,b,c,d) c./a;
persistence_y=@(a,b,c,d,base_a,base_b) (d-base_b)./(b-base_b);
persistence_x=@(a,b,c,d,base_a,base_b) (c-base_a)./(a-base_a);

useAssigns=unique(spikes.assigns);
grp_result=cell(1,length(grp_led));
for i=1:length(grp_led)
    curr_led=grp_led{i};
    curr_window=grp_window{i};
    curr_stimcond=grp_stimcond{i};
    grp_result{i}=getWindowFR(spikes,useAssigns,curr_stimcond,curr_led,curr_window,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED);
end

% ypoints=persistence_y(grp_result{4},grp_result{2},grp_result{3},grp_result{1});
% xpoints=persistence_x(grp_result{4},grp_result{2},grp_result{3},grp_result{1});
ypoints=persistence_y(grp_result{4},grp_result{2},grp_result{3},grp_result{1},grp_result{6},grp_result{5});
xpoints=persistence_x(grp_result{4},grp_result{2},grp_result{3},grp_result{1},grp_result{6},grp_result{5});

% ypoints(ypoints==0)=0.1;
% xpoints(xpoints==0)=0.1;
figure(); 
scatter(xpoints,ypoints);

sub=ypoints-xpoints;
[n,xout]=hist(sub(~isnan(sub) & ~isinf(sub)),20);
figure(); plot(xout,n);


end


function meanFRsForUnits=getWindowFR(spikes,useAssigns,useStimcond,useLedcond,window,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED)

meanFRsForUnits=zeros(length(useAssigns),1);
for i=1:length(useAssigns)
    useSpikes=filtspikes(spikes,0,'assigns',useAssigns(i),'stimcond',useStimcond);
    currTrialCon=unitByUnitTrials{i};
    currStimCon=unitByUnitStimcond{i};
    currLEDCon=unitByUnitLED{i};
    unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond) & ismember(currLEDCon,useLedcond));
    if isempty(unitByUnitConsensus)
        disp('unitByUnitConsensus is empty');
    end
    m=calcMeanAndStdForUnit(filtspikes(useSpikes,0,'led',useLedcond),window,unitByUnitConsensus);
    meanFRsForUnits(i)=m;
end
end

function [varargout]=calcMeanAndStdForUnit(spikes,window,theseTrials)

binsize=1; % in ms
if isempty(theseTrials)
    disp('theseTrials should not be empty in calcMeanAndStdForUnit');
end
numtrials=length(theseTrials);
% Set spiketimes
spiketimes = spikes.spiketimes;
% Convert binsize from ms to s
binsize = binsize/1000;
% Get counts
edges=window(1):binsize:window(2);
n = histc(spiketimes,edges);
n = n/numtrials/binsize;
m = mean(n);

nsForStdev=zeros(numtrials,size(n,2));
allTrials=theseTrials;
for i=1:length(allTrials)
    cspikes=filtpartialspikes(spikes,0,'trials',allTrials(i));
    if isempty(cspikes.spiketimes)
        continue
    end
    nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
end
nsForStdev=nsForStdev/binsize;
s=std(mean(nsForStdev,2),0,1);
if all(isnan(n))
    n = 0;
end
varargout{1} = mean(mean(nsForStdev,2));
varargout{2} = s;
varargout{3} = mean(nsForStdev,2);
end
    
