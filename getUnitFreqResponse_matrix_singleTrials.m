function [freqs,p2,responses,p,pSpont,bigFres]=getUnitFreqResponse_matrix_singleTrials(spikes,useAssign,leds)

usePowerSpecOfCrossCorr=1;
showFigs=0;
subtractSpont=0;
bigFres=[];
doChronux=1;

% freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% freqs=[1.03 2.03 4.03 6.03 8.03 10.03 12.03 14.03 16.03 18.03 20.03 30.03 40.03 50.03 60.03];
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% leds=[1.05 2.05 4.05 6.05 8.05 10.05 12.05 14.05 16.05 18.05 20.05 30.05 40.05 50.05 60.05];
% leds=[1.0 2.0 4.0 6.0 8.0 10.0 12.0 14.0 16.0 18.0 20.0 30.0 40.0 50.0 60.0];

p=zeros(length(leds),length(freqs));
pSpont=zeros(length(leds),length(freqs));
responses=zeros(15,513);
if showFigs==1
    figure(); 
end
for i=1:length(freqs)
%     disp(i);
    if doChronux==1
%         responses=[];
        p2=[];
        [~,~,~,~,unitsfre,unitsavSpec,~,~,unitsfreSpont,unitsavSpecSpont]=plot_unit_crossCorr_singleTrials(spikes,useAssign,leds(i),1);
        fre=unitsfre';
        avSpec=nanmean(unitsavSpec,2)';
        bigFres(i,:)=unitsfre(1,:);
        responses(i,:)=avSpec;
        for j=1:length(freqs)
            p(i,j)=max(avSpec(fre>=freqs(j)-0.5 & fre<=freqs(j)+0.5));
        end
    else
        for j=1:length(leds)
            if usePowerSpecOfCrossCorr==1
                % Uses power spectrum of cross-correlation of unit response with
                % stimulus
                %             tryS=makeTempField(spikes,'led',leds(j));
                %             temp(i,:)=tryS.temp;
                %             temp1(i,:)=tryS.sweeps.temp;
                %             tryS.temp=[];
                %             tryS.sweeps.temp=[];
                %             tryS.temp=sum(temp,1)>=1;
                %             tryS.sweeps.temp=sum(temp1,1)>=1;
                %             tryS=filtspikes(tryS,0,'temp',1,'assigns',useAssign);
                %             if isempty(tryS.led)
                %                 continue
                %             end
                [~,~,~,~,unitsfre,unitsavSpec,~,~,unitsfreSpont,unitsavSpecSpont]=plot_unit_crossCorr_singleTrials(spikes,useAssign,leds(j),freqs(i));
                %             if isempty(unitsfre)
                %                 continue
                %             end
                fre=mean(unitsfre,1);
                avSpec=mean(unitsavSpec,1);
                freSpont=mean(unitsfreSpont,1);
                avSpecSpont=mean(unitsavSpecSpont,1);
                if subtractSpont==1
                    if isempty(freSpont)
                        pSpont(j,i)=0;
                    else
                        pSpont(j,i)=max(avSpecSpont(freSpont>=1));
                    end
                end
            else
                % Uses power spectrum of auto-correlation of unit response
                [~,~,~,~,~,~,fre,avSpec,~,~,freSpont,avSpecSpont]=plot_unit_autoCorr(spikes,useAssign,leds(j));
                if subtractSpont==1
                    %avSpec=avSpec-avSpecSpont;
                    if isempty(freSpont)
                        pSpont(j,i)=0;
                    else
                        pSpont(j,i)=max(avSpecSpont(freSpont>=1));
                    end
                end
            end
            if isempty(fre) && i==j
                responses(i,:)=zeros(size(responses(1,:)));
                p(j,i)=0;
                continue
            end
            if i==j
                responses(i,:)=avSpec;
            end
            %         ma=max(avSpec(fre>=1));
            %         mi=min(avSpec(fre>=1));
            %         p(j,i)=max(avSpec(fre>=1));
            if isempty(fre)
                p(j,i)=0;
            else
                p(j,i)=max(avSpec(fre>=0.9));
            end
            if showFigs==0
                close all;
            else
                plot(fre,avSpec);
                hold on;
            end
        end
    end
end

if subtractSpont==1
    p2=p-pSpont;
else
    p2=p;
end

if showFigs==1
    figure();
    imagesc(responses);

    figure();
%     plot(freqs,p);
    imagesc(p);
    
end