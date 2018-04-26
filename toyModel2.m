function [time_steps,mua_output]=toyModel2()

% Generate spiking pdf for individual cells
% Exponential with cut-off decided by AP threshold
n_cells=100;
% time_steps=0:0.001:500; % in ms
% time_steps=0:0.1:2000; % in ms
time_steps=0:0.1:4000;
spiking_pdfs=zeros(n_cells,length(time_steps));
membrane_pdfs=zeros(n_cells,length(time_steps));
for i=1:n_cells
    if i>0 & i<=9 % L2/3
%         offset=13+0.5.*randn(1,1); % in ms
        offset=2.4+0.5.*randn(1,1); % in ms
        cell_tau=29; % in ms
        V_rest= -71; % in mV
        V_ap= -38.3; % in mV
        init_v = (0-V_rest) + 5.*randn(1,1); % Glu rev. at 0 mV (so mV above V_rest)
        thresh_v = (V_ap-V_rest) + 2.*randn(1,1); % mV above V_rest
    elseif i>9 & i<=46 % L4
%         offset=2+0.5.*randn(1,1); % in ms
        offset=2.4+0.5.*randn(1,1); % in ms
        cell_tau=34.8; % in ms
        V_rest= -66; % in mV
        V_ap= -39.7; % in mV
        init_v = (0-V_rest) + 5.*randn(1,1); % Glu rev. at 0 mV (so mV above V_rest)
        thresh_v = (V_ap-V_rest) + 2.*randn(1,1); % mV above V_rest
    elseif i>46 & i<=73 % L5a
%         offset=8+0.5.*randn(1,1); % in ms
        offset=2.4+0.5.*randn(1,1); % in ms
        cell_tau=37.6; % in ms
        V_rest= -62.8; % in mV
        V_ap= -38.9; % in mV
        init_v = (0-V_rest) + 5.*randn(1,1); % Glu rev. at 0 mV (so mV above V_rest)
        thresh_v = (V_ap-V_rest) + 2.*randn(1,1); % mV above V_rest
    elseif i>73 & i<=90 % L5b
%         offset=10+0.5.*randn(1,1); % in ms
        offset=2.4+0.5.*randn(1,1); % in ms
        cell_tau=25.8; % in ms
        V_rest= -63; % in mV
        V_ap= -41; % in mV
        init_v = (0-V_rest) + 5.*randn(1,1); % Glu rev. at 0 mV (so mV above V_rest)
        thresh_v = (V_ap-V_rest) + 2.*randn(1,1); % mV above V_rest
    elseif i>90 & i<=100 % L6
%         offset=10+0.5.*randn(1,1); % in ms
        offset=2.4+0.5.*randn(1,1); % in ms
        cell_tau=28.2; % in ms
        V_rest= -66.8; % in mV
        V_ap= -40.2; % in mV
        init_v = (0-V_rest) + 5.*randn(1,1); % Glu rev. at 0 mV (so mV above V_rest)
        thresh_v = (V_ap-V_rest) + 2.*randn(1,1); % mV above V_rest
    end
%     cell_tau=10;
    offset=floor(offset*10); % in indices
    exp_falloff=init_v.*exp(-time_steps./cell_tau);
    falloff_cut=[exp_falloff(1).*ones(1,offset) exp_falloff(1:end-offset)];
    falloff_cut(falloff_cut<thresh_v)=0;
    spiking_cut=exp_falloff-thresh_v;
    spiking_cut=[spiking_cut(1).*ones(1,offset) spiking_cut(1:end-offset)];
    spiking_cut(spiking_cut<0)=0;
    membrane_pdfs(i,:)=falloff_cut;
    spiking_pdfs(i,:)=(1/95)*spiking_cut;
    if i==10
        figure(); 
        plot(time_steps,membrane_pdfs(i,:));
        axis([0 100 0 100]);
        figure(); 
        plot(time_steps,spiking_pdfs(i,:));
        axis([0 100 0 1]);
    end
end

% Take average
figure(); 
plot(time_steps,mean(membrane_pdfs,1));
axis([0 100 0 100]);

figure(); 
mua_output=mean(spiking_pdfs,1);
plot(time_steps,mua_output);
hold all;
x=0:0.001:100;
plot(x,0.7.*exp(-x/10));
axis([0 100 0 0.7]);
