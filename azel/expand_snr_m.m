function [ snr_exp_m, snr_of_m ] = expand_snr_m(snr_m, varargin)
%EXPAND_SNR_M 
%   Inputs:
%       snr_m - Mx9 snr matrix
%       time_shift days to apply (optional)
%
%   Output:
%       snr_exp_m - expanded SNR matrix with one entry per second
%       snr_of_m  - overflow snr from time shift
%
%%%%%%%%%%%%%%%%%%
prn_max = 32;
snr_exp_m = zeros(prn_max,86400,9);

if nargin==2
    snr_m(:,1) = snr_m(:,1) - varargin{1}*snr_m(:,9);
    if varargin{1} < 0
        snr_of_m = snr_m(snr_m(:,1)>86400,:);
        snr_of_m(:,1) = snr_of_m(:,1) + 86400;
    elseif varargin{1} > 0
        snr_of_m = snr_m(snr_m(:,1)<=0,:);
        snr_of_m(:,1) = snr_of_m(:,1) + 86400;
    end
end

snr_m(:,1) = round(snr_m(:,1));
snr_m = snr_m(snr_m(:,1)>0,:);
snr_m = snr_m(snr_m(:,1)<=86400,:);

for jj = 1:prn_max
    snr_prn_m = snr_m(snr_m(:,2)==jj,:);
    snr_valid_idx = snr_prn_m(:,1)';
    snr_exp_m(jj,snr_valid_idx,:) = snr_prn_m;
end

end

