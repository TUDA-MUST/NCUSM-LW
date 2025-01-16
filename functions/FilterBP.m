function [bpFilter,bpFilter_f] = FilterBP(f_lo, f_up, N_filter, N, f_s)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    %% BP Filter#
 %Filter like this with convolution:
%  f_lo = 39000;
% f_up = 41000;
% N_filter = 20001;
% N = length(data);
% [bpFilter,bpFilter_f] = FilterBP(f_lo, f_up, N_filter, N, Fs);
% fil_data=conv(data, bpFilter);
% fil_data = fil_data(1:size(data,1),1:size(data,2));
% figure,title("check filtering"), hold on, plot(data./max(data)), plot(bpFilter./max(bpFilter)), plot(fil_data./max(fil_data))


n_up = ceil(N/f_s*f_up);
n_lo = floor(N/f_s*f_lo);

bpFilter_f = zeros(1,N);
bpFilter_f(1+n_lo:1+n_up) = 1;
bpFilter_f(end-n_up+1:end-n_lo+1) = 1;
bpFilter = ifft(bpFilter_f);
w_hann = hann(N_filter).';

bpFilter_trunc = zeros(1,N);
bpFilter_trunc(1:floor(N_filter/2)) = bpFilter(1:floor(N_filter/2)).*w_hann(round(N_filter/2+1):end); %watch out the floor and round does not really make sens for N odd
bpFilter_trunc(round(end-N_filter/2+1):end) = bpFilter(round(end-N_filter/2+1):end).*w_hann(1:floor(N_filter/2));

% bpFilter_causal = zeros(1,N);
% bpFilter_causal(1:N_filter/2) = bpFilter_trunc(end-N_filter/2+1:end);
% bpFilter_causal(N_filter/2+1:N_filter) = bpFilter_trunc(1:N_filter/2);

bpFilter = bpFilter_trunc;

bpFilter_f = (fft(bpFilter));
bpFilter= bpFilter./max(bpFilter_f);
bpFilter_f=bpFilter_f./max(bpFilter_f);

end

