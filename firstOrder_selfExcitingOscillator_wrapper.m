function firstOrder_selfExcitingOscillator_wrapper()


% tau=12; % in ms
% scale_rec=1;
% scale_inh=1.2;
% delta_t_step=4.12; 
% rec_delay=2;
% inh_delay=2;

% tau=12; % in ms
% scale_rec=1;
% scale_inh=1.2;
% delta_t_step=4.12; 
% rec_delay=2;
% inh_delay=2;

% tau=12; % in ms
% scale_rec=1;
% scale_inh=1;
% % delta_t_step=4.12; 
% delta_t_step=4;
% % rec_delay=3;
% % inh_delay=3;
% rec_delay=3.3;
% inh_delay=3;

tau=12; % in ms
scale_rec=1;
scale_inh=1;
% delta_t_step=4.12; 
delta_t_step=3;
% rec_delay=3;
% inh_delay=3;
rec_delay=3.4;
inh_delay=3;

try_tau=1:100;
try_ratioExctoInh=[0.01 0.05 0.1 0.25 0.5 1:50];
try_delta_t_step_asFractionOf_tau=[0.01:0.01:1 2:10];
try_diffBetweenExcAndInhDelay=[0.1:0.1:1 1:5];



firstOrder_selfExcitingOscillator(tau,scale_rec,scale_inh,delta_t_step,rec_delay,inh_delay);