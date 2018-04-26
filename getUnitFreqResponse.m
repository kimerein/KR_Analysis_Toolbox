function [freqs,p,responses]=getUnitFreqResponse(spikes,useAssign)

usePowerSpecOfCrossCorr=1;
showFigs=1;
subtractSpont=1;

% freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% freqs=[1.03 2.03 4.03 6.03 8.03 10.03 12.03 14.03 16.03 18.03 20.03 30.03 40.03 50.03 60.03];
freqs=[1.05 2.05 4.05 6.05 8.05 10.05 12.05 14.05 16.05 18.05 20.05 30.05 40.05 50.05 60.05];

p=zeros(1,length(freqs));
if showFigs==1
    figure(); 
end
for i=1:length(freqs)
    if usePowerSpecOfCrossCorr==1
        % Uses power spectrum of cross-correlation of unit response with
        % stimulus
        [~,~,~,~,fre,avSpec,~,~,freSpont,avSpecSpont]=plot_unit_autocorr(spikes,useAssign,freqs(i));
        if subtractSpont==1
            if isempty(freSpont)
                pSpont(i)=0;
            else
                pSpont(i)=max(avSpecSpont(freSpont>=1));
            end
        end
    else
        % Uses power spectrum of auto-correlation of unit response
        [~,~,~,~,~,~,fre,avSpec,~,~,freSpont,avSpecSpont]=plot_unit_autocorr(spikes,useAssign,freqs(i));
        if subtractSpont==1
            %avSpec=avSpec-avSpecSpont;
            if isempty(freSpont)
                pSpont(i)=0;
            else
                pSpont(i)=max(avSpecSpont(freSpont>=1));
            end
        end
    end
    if isempty(fre)
        responses(i,:)=zeros(size(responses(1,:)));
        p(i)=0;
        continue
    end
    responses(i,:)=avSpec;
    ma=max(avSpec(fre>=1));
    mi=min(avSpec(fre>=1));
    p(i)=max(avSpec(fre>=1));
    if showFigs==0
        close all;
    else
        plot(fre,avSpec);
        hold on;
    end
end

if subtractSpont==1
    p=p-pSpont;
end

if showFigs==1
    figure();
    imagesc(responses);

    figure();
    plot(freqs,p);
end