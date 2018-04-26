function firstOrder_selfExcitingOscillator(tau,scale_rec,scale_inh,delta_t_step,rec_delay,inh_delay)

% times=0:0.1:1000; % in ms
times=0:0.1:2000; % in ms
timestep=times(2)-times(1);

% tau=12; % in ms
% scale_rec=2;
% scale_inh=8;
% delta_t_step=2; 
% rec_delay=4;
% inh_delay=6;
% tau=12; % in ms
% scale_rec=2;
% scale_inh=50;
% delta_t_step=2; 
% rec_delay=4;
% inh_delay=6;
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
% rec_delay=2.8;
% inh_delay=3.2;
% tau=12; % in ms
% scale_rec=1;
% scale_inh=1.2;
% delta_t_step=4.12; 
% rec_delay=2.8;
% inh_delay=3;

% tau=12; % in ms
% scale_rec=1;
% scale_inh=2.11;
% delta_t_step=2; 
% rec_delay=8;
% inh_delay=10;

% tau=12; % in ms
% scale_rec=1;
% scale_inh=0.92;
% delta_t_step=2; 
% rec_delay=4;
% inh_delay=4.1;

% tau=12; % in ms
% scale_rec=1;
% scale_inh=100;
% delta_t_step=2; 
% rec_delay=1;
% inh_delay=1.1;
% scale_rec2=1;
% scale_inh2=0.7;
% rec_delay2=4;
% inh_delay2=4.1;

c1=1/tau;
x1_vals=zeros(1,length(times));
F_vals=zeros(1,length(times));
I_th=@(t) t>100 & t<900; % step function in ms
for i=200:length(times)
%     F=I_th(times(i))+scale_rec*max([x1_vals(find(times<times(i)-rec_delay*delta_t_step,1,'last')) 0])-scale_inh*max([x1_vals(find(times<times(i)-inh_delay*delta_t_step,1,'last')) 0]);
    F=I_th(times(i))+scale_rec*(x1_vals(find(times<times(i)-rec_delay*delta_t_step,1,'last')))-scale_inh*(x1_vals(find(times<times(i)-inh_delay*delta_t_step,1,'last')));
    F=scale_rec*(x1_vals(find(times<times(i)-rec_delay*delta_t_step,1,'last')))-scale_inh*(x1_vals(find(times<times(i)-inh_delay*delta_t_step,1,'last')));
%     F=I_th(times(i))+scale_rec*max([x1_vals(find(times<times(i)-rec_delay*delta_t_step,1,'last')) 0])-scale_inh*max([x1_vals(find(times<times(i)-inh_delay*delta_t_step,1,'last')) 0])+scale_rec2*max([x1_vals(find(times<times(i)-rec_delay2*delta_t_step,1,'last')) 0])-scale_inh2*max([x1_vals(find(times<times(i)-inh_delay2*delta_t_step,1,'last')) 0]);
    F_vals(i)=F;
%     x1_vals(i)=max([0 x1_vals(i-1)+timestep*(-c1*x1_vals(i-1)+F)]);
    x1_vals(i)=max([0 x1_vals(i-1)+I_th(times(i))+timestep*(-c1*x1_vals(i-1)+F)]);
%     x1_vals(i)=max([0 x1_vals(i-1)+timestep*(F)]);
end
upto=length(times);
figure();
subplot(3,1,1);
plot(times(1:upto),I_th(times(1:upto)));
xlabel('I_th');
subplot(3,1,2);
plot(times(1:upto),F_vals(1:upto));
xlabel('F_vals');
subplot(3,1,3);
plot(times(1:upto),x1_vals(1:upto));
xlabel('x1');