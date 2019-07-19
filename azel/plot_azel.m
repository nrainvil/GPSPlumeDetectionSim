function [ fig ] = plot_azel( data,az_bin,el_bin,data_title,ax_lim,varargin)

if(nargin<=5)
    ytxt = 'dB-Hz';
else
    ytxt = varargin{1};
end

%PLOT_AZEL
% Plot Azimuth/Elevation Plot
az_len = floor(360/az_bin);
el_len = floor(90/el_bin);
az_v = ((1:az_len+1)*az_bin);
el_v = (1:el_len)*el_bin;
el_v = [0,el_v];

%Translate data to pcolor format
data_plot = data;
data_plot = fliplr(data_plot')';
data_plot = fliplr(data_plot);
data_plot = ([data_plot,data_plot(:,end)]);
data_plot = fliplr(data_plot')';
data_plot = [data_plot;zeros(1,length(data_plot(1,:)))];
 
%Rotate, data is centered on AZ/EL spokes
az_v = az_v - az_bin/2;
fig = figure;
[h,c] = polarPcolor(el_v,az_v,data_plot,...
            'Ncircles',7,'Nspokes',9);
title(data_title);
caxis(ax_lim);
ylabel(c,ytxt);

end