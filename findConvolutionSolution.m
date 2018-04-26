function [f_x,f_y,f_v,get_alternate]=findConvolutionSolution(f_x,f_y,f_v)
% f_v=convolution of f_x and f_y
% if f_v is missing, find f_v
% if f_x is missing, find f_x
% if f_y is missing, find f_y

% Example
% r1 = 75 + 50.*randn(1000,1);
% r2 = 80 + 40.*randn(1000,1);
% [n,x]=hist(r1+r2,100);
% pdf_r1=hist(r1,x);
% pdf_r2=hist(r2,x);
% c=conv(pdf_r1,pdf_r2,'same');
% figure(); plot(x,n*1000); hold all; plot(x+130,c);
% [d,r]=deconv(c',pdf_r1');
% figure(); plot(x,(d-r)./10000000); hold all; plot(x,pdf_r2);

if isempty(f_v)
    % Just get convolution of f_x and f_y
    N=length(f_x)+length(f_y)-1;
    fft_fx=fft(f_x,N);
    fft_fy=fft(f_y,N);
    fft_fv=fft_fx.*fft_fy;
    f_v=real(ifft(fft_fv));
    get_alternate=conv(f_x,f_y);
end
if isempty(f_x)
    % Deconvolve
    N=length(f_v)+length(f_y)-1;
    fft_fv=fft(f_v,N);
    fft_fy=fft(f_y,N);
    fft_fx=fft_fv./fft_fy;
    f_x=real(ifft(fft_fx));
    get_alternate=deconv(f_v,f_y);
%     [q,r]=deconv(f_v,f_y);
%     f_x=q;
%     get_alternate=r;
end
if isempty(f_y)
    % Deconvolve
    [q,r]=deconv(f_v,f_x);
    f_y=q;
    get_alternate=r;
end



