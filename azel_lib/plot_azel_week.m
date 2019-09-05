function [fig1, fig2] = plot_azel_week( rcvr, start_day, year, save_l)
%PLOT_AZEL_WEEK 
%Plot AZ EL Statistics for a week
az_bin = 45/2;
el_bin = 5;

end_day = start_day+6;
day_v = mod(start_day:end_day,365);
day_v(day_v==0) = 365;
year_v = year.*ones(1,length(day_v));
year_v(day_v<start_day) = year+1;

[snr_all_m,az_all_m,el_all_m,ts_all_m] = load_snr_all_m(rcvr,day_v,year_v);

az_len = floor(360/az_bin);
el_len = floor(90/el_bin);
snr_azel = zeros(az_len,el_len);
for ii = 1:az_len
    for jj = 1:el_len
        az = ii*az_bin;
        el = jj*el_bin;
        %fprintf('AZ: %03.2f (%d) EL: %02.2f (%d) \n',az,el,ii,jj);
        ii_shift = mod(ii + 90/az_bin,az_len);
        
        snr_ele = snr_all_m(bin_deg(az_all_m,az,az_bin));
        el_ele =  el_all_m(bin_deg(az_all_m,az,az_bin));
        snr_ele_polar = snr_ele(bin_deg(el_ele,el,el_bin));
        snr_azel_polar(ii,jj) = mean(snr_ele_polar);
        std_azel_polar(ii,jj) = std(snr_ele_polar);
    end
end

fig1 = plot_azel(snr_azel_polar,az_bin,el_bin,...
                    strcat(upper(rcvr),' Mean SNR'),[20 55]);
fig2 = plot_azel(std_azel_polar,az_bin,el_bin,...
                    strcat(upper(rcvr),' Std-Dev SNR'),[0 5]);

if save_l
    fname = strcat(rcvr,'_',num2str(year),...
                        '_',num2str(start_day,'%03.0f'),...
                        '_',num2str(mod(end_day,365),'%03.0f'));
    print(fig1,strcat('auto_plots\snr_azel_',fname),'-dpng');
    print(fig2,strcat('auto_plots\std_azel_',fname),'-dpng');
end

end

