function toyConductanceModel()

cell_tau=12; % in ms
thresh_v_range=25:0.001:45; % -40 mV to -20 mV, use -65 as V_rest
init_range=35:0.001:95; % max. of dendritic depolarization
init_pdf=@(x) normpdf(x,65,10); % Glu rev. at 0 mV
thresh_pdf=@(x) normpdf(x,35,5); % Mean at -30 mV, variance = 5 mV?
% Let E and I conductances convolved with cell_tau direct membrane potential
tau_E=1.5; % in ms
tau_I=5; % in ms
EI_ratio=5;

% figure(); 
% plot(thresh_v_range,thresh_pdf(thresh_v_range));

% Generate spiking pdf for individual cells
% Exponential with cut-off decided by AP threshold
n_cells=5;
time_steps=0:0.001:100; % in ms
spiking_pdfs=zeros(n_cells,length(time_steps));
membrane_pdfs=zeros(n_cells,length(time_steps));
for i=1:n_cells
    disp(i);
    init_v = 65 + 30.*randn(1,1);
    thresh_v = 20 + 12.*randn(1,1);  
    ei_r = EI_ratio + 0.1.*randn(1,1);
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
    offset=floor(offset*1000); % in indices
%     exp_falloff=init_v.*exp(-time_steps./cell_tau);
    % Net EI effects
    add_inh=exp(-time_steps./tau_I);
    add_inh(1)=0;
    net=ei_r.*exp(-time_steps./tau_E)-add_inh;
    m_cap=exp(-time_steps./cell_tau);
    m_timeCourse=conv(net,m_cap,'same');
    scaleFactor=init_v./m_timeCourse(1);
    m_timeCourse=m_timeCourse.*scaleFactor;
    exp_falloff=m_timeCourse;
    falloff_cut=[exp_falloff(1).*ones(1,offset) exp_falloff(1:end-offset)];
    firstoff=find(falloff_cut<thresh_v,1,'first');
    falloff_cut(firstoff:end)=0;
    spiking_cut=exp_falloff-thresh_v;
    spiking_cut=[spiking_cut(1).*ones(1,offset) spiking_cut(1:end-offset)];
    firstoff=find(spiking_cut<0,1,'first');
    spiking_cut(firstoff:end)=0;
    membrane_pdfs(i,:)=falloff_cut;
    spiking_pdfs(i,:)=(1/95)*spiking_cut;
    if i==1
        figure(); 
        plot(time_steps,net);
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
plot(time_steps,mean(spiking_pdfs,1));
hold all;
x=0:0.001:100;
plot(x,0.7.*exp(-x/10));
