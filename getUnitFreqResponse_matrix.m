function [freqs,p2,responses,p,pSpont,bigFres]=getUnitFreqResponse_matrix(spikes,useAssign,leds)

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
if showFigs==1
    figure(); 
end
responses=[];
bigFres=[];
p2=[];
for i=1:length(freqs)
%     disp(i);
    if doChronux==1
        [~,~,~,~,fre,avSpec,~,~,freSpont,avSpecSpont]=plot_unit_crossCorr(spikes,useAssign,leds(i),1);
        bigFres(i,:)=fre;
        responses(i,:)=avSpec;
        for j=1:length(freqs)
%             if i==12 & j==12
%                 disp('hi');
%             end
            p(i,j)=max(avSpec(fre>=freqs(j)-0.5 & fre<=freqs(j)+0.5));
%             p(i,j)=max(avSpec(fre>=freqs(j)-1 & fre<=freqs(j)+1));
        end
    else
        for j=1:length(leds)
            %         disp(j);
            if usePowerSpecOfCrossCorr==1
                % Uses power spectrum of cross-correlation of unit response with
                % stimulus
                [~,~,~,~,fre,avSpec,~,~,freSpont,avSpecSpont]=plot_unit_crossCorr(spikes,useAssign,leds(j),freqs(i));
                if subtractSpont==1
                    if isempty(freSpont)
                        pSpont(j,i)=0;
                    else
                        pSpont(j,i)=max(avSpecSpont(freSpont>=1));
                    end
                end
            else
                useAssign=[];
                % Uses power spectrum of auto-correlation of unit response
                [~,~,~,~,~,~,fre,avSpec,~,~,freSpont,avSpecSpont]=plot_unit_autocorr(spikes,useAssign,leds(j));
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
            if usePowerSpecOfCrossCorr==1 && i==j
                responses(i,:)=avSpec;
            elseif usePowerSpecOfCrossCorr==0
                responses(j,:)=avSpec;
                bigFres=fre;
            end
            %         ma=max(avSpec(fre>=1));
            %         mi=min(avSpec(fre>=1));
            p(j,i)=max(avSpec(fre>=1));
            p(j,i)=max(avSpec(fre>=freqs(i)-0.5 & fre<=freqs(i)+0.5));
            %         p(j,i)=max(avSpec(fre>=0.9));
            %         p(j,i)=max(avSpec(fre>=0.97));
            if showFigs==0
                close all;
            else
                plot(fre,avSpec);
                hold on;
            end
        end
        if usePowerSpecOfCrossCorr==0 && j==length(leds)
            break
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