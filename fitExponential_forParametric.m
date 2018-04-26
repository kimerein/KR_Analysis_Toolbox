function [finalTau,finalW,finalResidNorm,f]=fitExponential_forParametric(x,y,LEDonset,fitTau,waitS,tauDefault,zeroWindow,showTestFigs)
% x and y are row vectors
% LEDonset is the time in seconds when the LED comes on

% Fit to the following function form:
% y=A*exp(-t/tau)

alreadyNormalized=1;
w_range=0.000:0.001:0.03;
% zeroWindow=[0.05 0.08];

if fitTau==0
    % Fit w rather than tau
    % Use default tau
    residNorms=zeros(length(w_range),1);
    for i=1:length(w_range)
        waitS=w_range(i);
%         finishS=waitS+zeroWindow(2);
        finishS=zeroWindow(2);
        subX=x(x>=LEDonset+waitS & x<=LEDonset+finishS);
        subX=subX-min(subX);
        subY=y(x>=LEDonset+waitS & x<=LEDonset+finishS);
        if alreadyNormalized==0
            % Zero curve during zeroWindow
            subY=subY-mean(subY(subX>=zeroWindow(1) & subX<=zeroWindow(2)));
            % Normalize top of curve to 1
            topVal=mean(subY(1:find(subX==LEDonset)));
            subY=subY./topVal;
        end
        % Get log of subY
        % y-intercept must be 0, because normalized
        % so A=1
        % slope of linear fit will give tau
        resid=subY-exp(-subX/tauDefault);
        residNorms(i)=(sum(resid.^2))^(1/2);
    end
%     disp([residNorms w_range']);
    [mi,minind]=min(residNorms,[],1);
    finalTau=tauDefault;
    finalW=w_range(minind);
    finalResidNorm=mi;
    subX=x(x>=LEDonset+finalW & x<=LEDonset+finishS);
    subX=subX-min(subX);
    subY=y(x>=LEDonset+finalW & x<=LEDonset+finishS);
else
    % Fit tau rather than w
    % Use waitS as default
     finishS=waitS+zeroWindow(2);
     subX=x(x>=LEDonset+waitS & x<=LEDonset+finishS);
     subX=subX-min(subX);
     subY=y(x>=LEDonset+waitS & x<=LEDonset+finishS);
     if alreadyNormalized==0
         % Zero curve during zeroWindow
         subY=subY-mean(subY(subX>=zeroWindow(1) & subX<=zeroWindow(2)));
         % Normalize top of curve to 1
         topVal=mean(subY(1:find(subX==LEDonset)));
         subY=subY./topVal;
     end
     firstZero=find(subY<=0,1,'first');
     % Get log of subY
     log_subY=log(subY(1:firstZero-1)); % linear if y is exponential
     log_subX=subX(1:firstZero-1);
     % Fit a line to log_subY
     p=polyfit(log_subX',log_subY',1);
     % Get exponential fit from linear fit
     % p(1) is the slope, p(2) is the y-intercept
     tau=-1/p(1);
     A=exp(p(2));
     resid=subY-A*exp(-subX/tau);
     finalTau=tau;
     finalW=waitS;
     finalResidNorm=(sum(resid.^2))^(1/2);
end
 

% Plot results
f=[];
if showTestFigs==1
    f=figure();
    plot(subX,subY,'Color','b');
    hold on;
    plot(subX,exp(-subX/finalTau),'Color','r');
end
