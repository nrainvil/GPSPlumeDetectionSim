function [ snr_all_m, az_all_m, el_all_m,ts_all_m ] = ...
                                load_snr_all_m(rcvr,doy_v,year_v,varargin)
%LOAD_SNR_ALL_M 
prn_max=32;
time_shift = 0;
freq = 1;   %L1=1,L2=2
print_loading = 0;
force_load = 0;
if nargin>=4
    freq = varargin{1};
end
if nargin>=5
    force_load = varargin{2};
end

%Check Inputs
validateattributes(rcvr,{'char','cell'},{'nonempty'});
validateattributes(doy_v,{'numeric'},{'nrows',1});
day_cnt = length(doy_v(1,:));
validateattributes(year_v,{'numeric'},{'nrows',1,'ncols',day_cnt});

snr_all_m = zeros(day_cnt,prn_max,86400);
az_all_m = zeros(day_cnt,prn_max,86400);
el_all_m = zeros(day_cnt,prn_max,86400);
ts_all_m = zeros(day_cnt,prn_max,86400);

%Load Each Day
for day_idx=1:length(doy_v(1,:))
    doy = doy_v(day_idx);
    year = year_v(day_idx);
    if(doy > 365)
        doy = mod(doy,365);
        year = year+1;
    end
    snr_m = load_single_snr(rcvr,doy,year,print_loading,force_load);
    if(time_shift)
        snr_exp_m = expand_snr_m(snr_m,length(doy_v(1,:))-day_idx);
    else
        snr_exp_m = expand_snr_m(snr_m);
    end
    el_all_m(day_idx,:,:) = squeeze(snr_exp_m(:,:,3));
    az_all_m(day_idx,:,:) = squeeze(snr_exp_m(:,:,4));
    snr_all_m(day_idx,:,:) = squeeze(snr_exp_m(:,:,4+freq));
    ts_all_m(day_idx,:,:) = squeeze(snr_exp_m(:,:,9));
end

end

