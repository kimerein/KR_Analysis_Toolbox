function [all_thalSet,all_V1Set,concat_thalSet,concat_V1Set]=getTransferFunction(xpoints,thal,V1,freqs,diffStimOnsets,concat)

sOnset=1;
cycles=1:40;
if diffStimOnsets==1
    for i=1:length(freqs)
        stimOnsets{i}=(1/freqs(i)).*(cycles-1)+sOnset;
    end
else
    stimOnsets=1;
end

thalDelay=0.075;
V1Delay=0.075;
decimateRate=1;
smoothBin=2/2;
colors={[0.1 0.1 0.1]; [0.8 0.1 0.1]; [0.1 0.8 0.1]; [0.1 0.1 0.8]; ...
        [0.5 0 0.5]; [0 0.5 0.5]; [0.8 0.8 0.8]};

all_thalSet=cell(length(freqs),1);
all_V1Set=cell(length(freqs),1);
concat_thalSet=[];
concat_V1Set=[];
figure(); 
for i=1:length(freqs)
    if diffStimOnsets==0
        fourthPeriod=((1/freqs(i))/4)*4;
        thalSet1=smooth(thal(i,xpoints>stimOnsets+thalDelay & xpoints<stimOnsets+thalDelay+fourthPeriod),smoothBin);
        thalSet=thalSet1(1:decimateRate:end);
        V1Set1=smooth(V1(i,xpoints>stimOnsets+V1Delay & xpoints<stimOnsets+V1Delay+fourthPeriod),smoothBin);
        V1Set=V1Set1(1:decimateRate:end);
        m=min(length(thalSet),length(V1Set));
        plot(thalSet(1:m),V1Set(1:m),'Color',colors{i});
        line([0 3],[300-i*2 300-i*2],'Color',colors{i});
        hold on;
        all_thalSet{i}=thalSet;
        all_V1Set{i}=V1Set;
        concat_thalSet=[concat_thalSet; thalSet];
        concat_V1Set=[concat_V1Set; V1Set];
    else
        theseCycles=stimOnsets{i};
        for j=1:length(theseCycles)
            fourthPeriod=((1/freqs(i))/4);
            thalSet1=smooth(thal(i,xpoints>theseCycles(j)+thalDelay & xpoints<theseCycles(j)+thalDelay+fourthPeriod),smoothBin);
            thalSet=thalSet1(1:decimateRate:end);
            V1Set1=smooth(V1(i,xpoints>theseCycles(j)+V1Delay & xpoints<theseCycles(j)+V1Delay+fourthPeriod),smoothBin);
            V1Set=V1Set1(1:decimateRate:end);
            m=min(length(thalSet),length(V1Set));
            thalSet=thalSet(1:m);
            V1Set=V1Set(1:m);
            if j==1
                av_thalSet=thalSet;
                av_V1Set=V1Set;
            else
                if length(thalSet)==0
                    break
                end
                if concat==0
                    if length(thalSet)>length(av_thalSet)
                        thalSet=thalSet(1:length(av_thalSet));
                        V1Set=V1Set(1:length(av_V1Set));
                    elseif length(thalSet)<length(av_thalSet)
                        av_thalSet=av_thalSet(1:length(thalSet));
                        av_V1Set=av_V1Set(1:length(V1Set));
                    end
                    av_thalSet=av_thalSet+thalSet;
                    av_V1Set=av_V1Set+V1Set;
                else
                    av_thalSet=[av_thalSet; thalSet];
                    av_V1Set=[av_V1Set; V1Set];
                end
            end
        end
        if concat==0
            av_thalSet=av_thalSet/j;
            av_V1Set=av_V1Set/j;
        end
        plot(av_thalSet,av_V1Set,'Color',colors{i});
        line([0 3],[300-i*2 300-i*2],'Color',colors{i});
        hold on;
        all_thalSet{i}=av_thalSet;
        all_V1Set{i}=av_V1Set;
        concat_thalSet=[concat_thalSet; av_thalSet];
        concat_V1Set=[concat_V1Set; av_V1Set];
    end
end
    
