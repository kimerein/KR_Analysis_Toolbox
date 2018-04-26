function [fraction_under,fraction_over]=get_under_and_over_estimates(output)

% Under-estimate
% 1. p-val less than 0.05
% 2. amplitude change above 0
% 3. suppression of evoked response above 0
% 4. HWHM above 0.22 * 10^-3
% 5. include units at all depths
% 
% Over-estimate
% 1. p-val less than 0.05 OR 
% p-val less than or equal to 1 AND suppression of evoked response above 0
% 2. amplitude change above 0
% 3. HWHM above 0.22 * 10^-3
% 4. include only units below first putative Ntsr1+ unit

% output = 
% 
%       amp_change: [1x178 double]
%            pvals: [1x178 double]
%     visev_change: [1x178 double]
%           colors: [178x3 double]
%           depths: [178x1 double]
%       halfWidths: [178x1 double]
%            wvfms: [178x38 double]

output.depths=output.depths';
output.halfWidths=output.halfWidths';

% Under-estimate
true1=output.pvals<0.05;
true2=output.amp_change>0;
true3=output.visev_change>0;
true4=output.halfWidths>0.22*10^-3;
n_under=sum(true1 & true2 & true3 & true4);
includeUnits_under=sum(ones(size(output.depths)) & output.halfWidths>0.22*10^-3);
fraction_under=n_under/includeUnits_under;

% Over-estimate
true1=output.pvals<0.05;
true1Alternate=output.pvals<=1;
true2=output.amp_change>0;
true3=output.visev_change>0;
true4=output.halfWidths>0.22*10^-3;
temp=(true1 | (true1Alternate & true3)) & true2 & true4;
n_over=sum(temp);
includeUnits_over=sum(output.depths>min(output.depths(temp)) & output.halfWidths>0.22*10^-3);
fraction_over=n_over/includeUnits_over;
