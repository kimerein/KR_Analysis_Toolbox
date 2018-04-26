function selfExcitingOscillator()

% dampingRatio=0.5; 
% % w_0=2*pi*(60/1000); 
% tau=0.3*4.12; % in ms
% w_0=1/tau;
% delta_t_step=4.12; 
% scale_rec=0.25;
% scale_inh=0.6;

% dampingRatio=2.8; 
% w_0=2*pi*(50/1000); 
% % tau=0.3*4.12; % in ms
% % w_0=1/tau;
% delta_t_step=4.12; 
% scale_rec=1;
% scale_inh=1;

% dampingRatio=1.1; 
% w_0=2*pi*(150/1000); 
% % tau=0.3*4.12; % in ms
% % w_0=1/tau;
% delta_t_step=4.12; 
% scale_rec=1;
% scale_inh=1.5;

% dampingRatio=1; 
% w_0=2*pi*(150/1000); 
% % tau=0.3*4.12; % in ms
% % w_0=1/tau;
% delta_t_step=4.12; 
% scale_rec=1;
% scale_inh=2;

% dampingRatio=1; 
% w_0=2*pi*(150/1000); 
% % tau=0.3*4.12; % in ms
% % w_0=1/tau;
% delta_t_step=4.12; 
% scale_rec=1;
% scale_inh=0.1;

% dampingRatio=10; 
% w_0=2*pi*(30/1000); 
% % tau=0.3*4.12; % in ms
% % w_0=1/tau;
% delta_t_step=8; 
% scale_rec=2;
% scale_inh=2.2;

% dampingRatio=10; 
% w_0=2*pi*(30/1000); 
% % tau=0.3*4.12; % in ms
% % w_0=1/tau;
% delta_t_step=2; 
% scale_rec=1;
% scale_inh=2.11;
% rec_delay=8;
% inh_delay=10;

dampingRatio=10; 
w_0=2*pi*(150/1000); 
% tau=0.3*4.12; % in ms
% w_0=1/tau;
delta_t_step=4.12; 
scale_rec=4;
scale_inh=5;
rec_delay=2;
inh_delay=1.6;

times=0:0.1:1000; % in ms
timestep=times(2)-times(1);
v_vals=zeros(1,length(times));
x_vals=zeros(1,length(times));
F_vals=zeros(1,length(times));
I_rec_vals=zeros(1,length(times));
I_inh_vals=zeros(1,length(times));
I_th=@(t) t>100 & t<900; % step function input from thalamus
for i=200:length(times)
    F_vals(i)=I_th(times(i))+scale_rec*max([x_vals(find(times<times(i)-rec_delay*delta_t_step,1,'last')) 0])-scale_inh*max([x_vals(find(times<times(i)-inh_delay*delta_t_step,1,'last')) 0]);
%     F_vals(i)=I_th(times(i))+scale_rec*max([x_vals(find(times<times(i)-1*delta_t_step,1,'last')) 0])+scale_rec*max([x_vals(find(times<times(i)-2*delta_t_step,1,'last')) 0])+scale_rec*max([x_vals(find(times<times(i)-3*delta_t_step,1,'last')) 0])-scale_inh*max([x_vals(find(times<times(i)-4*delta_t_step,1,'last')) 0]);
    v_vals(i)=v_vals(i-1)+timestep*(F_vals(i)-w_0^2*x_vals(i-1)-2*dampingRatio*w_0*v_vals(i-1));
    x_vals(i)=max([0 x_vals(i-1)+timestep*v_vals(i-1)]);
    I_rec_vals(i)=scale_rec*x_vals(find(times<times(i)-2*delta_t_step,1,'last'));
    I_inh_vals(i)=-scale_inh*x_vals(find(times<times(i)-1.6*delta_t_step,1,'last'));
end
upto=length(times);
figure();
subplot(6,1,1);
plot(times(1:upto),I_th(times(1:upto)));
xlabel('I_th');
subplot(6,1,2);
plot(times(1:upto),F_vals(1:upto));
xlabel('F');
subplot(6,1,3);
plot(times(1:upto),I_rec_vals(1:upto));
xlabel('I_rec');
subplot(6,1,4);
plot(times(1:upto),I_inh_vals(1:upto));
xlabel('I_inh');
subplot(6,1,5);
plot(times(1:upto),x_vals(1:upto));
xlabel('x');
subplot(6,1,6);
plot(times(1:upto),v_vals(1:upto));
xlabel('v');







% solvedX=ones(1,length(times)).*x_init;
% dx_over_dt=zeros(1,length(times));
% keep_d2x_over_dt2=zeros(1,length(times));
% 
% x = @(t) solvedX(find(times<t,1,'last'));
% I_th= @(t) t>100 & t<900; % step function input from thalamus
% I_rec = @(t) x(t-2*delta_t_step);
% I_inh = @(t) -x(t-1.6*delta_t_step);
% 
% for i=200:length(times)
%     if i==1
%         d2x_over_dt2=0;
%     else
%         d2x_over_dt2=(dx_over_dt(i-1)-dx_over_dt(i-2))./(times(i-1)-times(i-2));
%         keep_d2x_over_dt2(i)=d2x_over_dt2;
%     end
%     F=I_th(times(i-1))+scale_rec*I_rec(times(i-1))+scale_inh*I_inh(times(i-1));
%     delta_x=(F-w_0^2*x(times(i-1))-d2x_over_dt2)*(1/(2*dampingRatio*w_0))*(times(i)-times(i-1));
%     solvedX(i)=x(times(i-1))+delta_x;
%     dx_over_dt(i)=delta_x/(times(i)-times(i-1));
%     if isinf(delta_x)
%         disp('hi');
%     end
% end
% 
% figure(); 
% subplot(3,1,1);
% plot(times,solvedX);
% subplot(3,1,2);
% plot(times,dx_over_dt);
% subplot(3,1,3);
% plot(times,keep_d2x_over_dt2);
%     
%         