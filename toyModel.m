function [time_steps,mua_output]=toyModel()

cell_tau=12; % in ms
thresh_v_range=25:0.001:45; % -40 mV to -20 mV, use -65 as V_rest
init_range=35:0.001:95; % max. of dendritic depolarization
init_pdf=@(x) normpdf(x,65,10); % Glu rev. at 0 mV
thresh_pdf=@(x) normpdf(x,35,5); % Mean at -30 mV, variance = 5 mV?

% figure(); 
% plot(thresh_v_range,thresh_pdf(thresh_v_range));

% Generate spiking pdf for individual cells
% Exponential with cut-off decided by AP threshold
n_cells=100;
% time_steps=0:0.001:500; % in ms
% time_steps=0:0.1:2000; % in ms
time_steps=0:0.1:4000;
spiking_pdfs=zeros(n_cells,length(time_steps));
membrane_pdfs=zeros(n_cells,length(time_steps));
for i=1:n_cells
    init_v = 65 + 30.*randn(1,1);
    thresh_v = 20 + 12.*randn(1,1);  
    offset=-1;
    while offset<0
%         offset=17+10.*randn(1,1); % in ms
%         offset=5+3.*randn(1,1); % in ms
        if i>0 & i<=9 % L2/3
            offset=13+0.5.*randn(1,1); % in ms
        elseif i>9 & i<=46 % L4
            offset=2+0.5.*randn(1,1); % in ms
        elseif i>46 & i<=73 % L5a
            offset=8+0.5.*randn(1,1); % in ms   
        elseif i>73 & i<=90 % L5b
            offset=10+0.5.*randn(1,1); % in ms 
        elseif i>90 & i<=100 % L6
            offset=10+0.5.*randn(1,1); % in ms 
        end
    end
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
        figure(); 
        plot(time_steps,spiking_pdfs(i,:));
    end
end

% Take average
figure(); 
plot(time_steps,mean(membrane_pdfs,1));
figure(); 
mua_output=mean(spiking_pdfs,1);
plot(time_steps,mua_output);
hold all;
x=0:0.001:100;
plot(x,0.7.*exp(-x/10));
