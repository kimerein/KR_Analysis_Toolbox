function toyTFtuning()

% If every neuron is a first-order low-pass filter with a different 
% cut-off frequency ... what does MUA look like with respect to single-unit
% TF tuning?

n=100; % Try 100 cells
x=0.5:0.001:40;
av_y=zeros(n,length(x));
for i=1:n
    curr_f_c=-1;
    while curr_f_c<0
        curr_f_c=9+3.*randn(1,1);
    end
    curr_y=amp_bodeplot(curr_f_c,x);
    av_y(i,:)=curr_y;
end

figure(); 
plot(x,mean(av_y,1));

end

function y=amp_bodeplot(f_c,x)

y=zeros(1,length(x));
for i=1:length(x)
    if x(i)<=f_c
        y(i)=1; % No reduction in amplitude
    else
        % How many times has frequency doubled
        n=log2(x(i)/f_c);
        % Gives how many times amplitude reduced by half
        y(i)=1/(2^n);
    end
end
end