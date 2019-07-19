% close all;
clear all;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('lib');
addpath('emag_lib');
addpath('plotGen');

%Setup
plume_rad = 500;
plume_height = 5000; %m

%Receivers
ecne = VaprRcvr(37+45/60+55.08/3600,15+6.48/3600,2892,'ECNE');
eplu = VaprRcvr(37+45/60+54.36/3600,14+59/60+8.52/3600,2913,'EPLU');
epdn = VaprRcvr(37+45/60+57.24/3600,15+1/60+.48/3600,2804,'EPDN');
esln = VaprRcvr(37+41/60+36/3600,14+58/60+27/3600,1728,'ESLN');
ecpn = VaprRcvr(37+44/60+37/3600,14+59/60+11/3600,2991,'ECPN');

%Vents
nec = Location(37+45/60+17.19/3600,14+59/60+48.41/3600,3292);
vor = Location(37+45/60+2.33/3600,14+59/60+36.82/3600,3267);
sec = Location(37+44/60+50.17/3600,15+3.67/3600,3257);

center = nec;
vents_c = {nec,vor,sec};
vents_txt_c = {'NEC','VOR','SEC'};
% vents_c = {sec};
% vents_txt_c = {'SEC'};
rcvr_c = {ecne,eplu,epdn,esln,ecpn};
rcvr_txt_c = {'RX13','RX14','RX16','RX11','RX15'};
% rcvr_c = {ecne,ecpn};
% rcvr_txt_c = {'RX13','RX15'};

figure;
hold on;
set(gca,'Color',[.9,.9,.9])
for rIdx=1:length(rcvr_txt_c)
    rcvr = rcvr_c{rIdx};
    [rENU(1), rENU(2), rENU(3)] = rcvr.enu(center);
    el_vent = [];
for vIdx=1:length(vents_c)
    %Plot Vents
    vent = vents_c{vIdx};
    [vENU(1,vIdx), vENU(2,vIdx), vENU(3,vIdx)] = vent.enu(center);
    [az(:,vIdx),vTang{vIdx}] = findAzLim( vENU(:,vIdx)', rENU, plume_rad);
    vent_az_c{vIdx} = 4000*[cosd(az(1,vIdx)),sind(az(1,vIdx)); ...
                           cosd(az(2,vIdx)),sind(az(2,vIdx))];
    %Plot AZ Limits
    viscircles([vENU(1,vIdx),vENU(2,vIdx)],plume_rad,'LineStyle',':');

    %Plot Vents
    vent = vents_c{vIdx};
    [vENU(1,vIdx), vENU(2,vIdx), vENU(3,vIdx)] = vent.enu(center);
    plot(vENU(1,vIdx),vENU(2,vIdx),'dr','LineWidth',2);
    text(vENU(1,vIdx)+50,vENU(2,vIdx)+75,vents_txt_c{vIdx},'Color','r');

    azel = [200,14];
    [edge_dist, delta_height, ppENU ] = ...
        findPP( vent, rcvr, azel, plume_rad, plume_height, center);
    if(edge_dist>0)
        az_adj = matAz(azel(1));
        az_ref_v = edge_dist*[cosd(az_adj),sind(az_adj)]; %Direction vector
        fprintf('%s intersects %s,',rcvr.name,vents_txt_c{vIdx});
        fprintf(' Dist: %0.2f m Height: %0.2f m\n',...
                edge_dist,delta_height);
        line([rENU(1),rENU(1)+az_ref_v(1)],...
             [rENU(2),rENU(2)+az_ref_v(2)],...
             'LineWidth',3,'LineStyle','-','Color',[0 1 0]);
    end
end

az_c{rIdx} = az;
end

for rIdx=1:length(rcvr_txt_c)
    rcvr = rcvr_c{rIdx};
    [rENU(1), rENU(2), rENU(3)] = rcvr.enu(center);
%Plot Max AZ Limit
az = az_c{rIdx};
az_max = max(max(az));
az_min = min(min(az));
az_max_v = 4000*[cosd(az_max),sind(az_max)]; 
az_min_v = 4000*[cosd(az_min),sind(az_min)]; 
line([rENU(1),rENU(1)+az_max_v(1)],[rENU(2),rENU(2)+az_max_v(2)],...
        'LineWidth',2);
line([rENU(1),rENU(1)+az_min_v(1)],[rENU(2),rENU(2)+az_min_v(2)],...
    'LineWidth',2);


fprintf('%s FOV: %0.2f to %0.2f\n',rcvr_txt_c{rIdx},gpsAz(az_max),gpsAz(az_min));

%Plot Receivers
plot(rENU(1),rENU(2),'ok','LineWidth',2,'MarkerSize',10);
text(rENU(1)+75,rENU(2)+75,rcvr_txt_c{rIdx});
end

xlim([-1500,1500]);
ylim([-1500,1500]);
set(gcf,'Position',[100 100 500 500])
title('Field of view for Etna GPS Stations');
xlabel('East (m)');
ylabel('North (m)');
xticks([-2000,-1000,0,1000,2000]);
yticks([-2000,-1000,0,1000,2000]);