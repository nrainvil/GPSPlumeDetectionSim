clear all;
% close all;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('lib');
addpath('emag_lib');
addpath('plotGen');

minEle = 20; %Degrees
antHeight = 2; % m 

%Simulation time 
startTime = datetime(2013,11,23,0,0,0,'TimeZone','UTC');
endTime = datetime(2013,11,23,23,55,0,'TimeZone','UTC');
deltaTime = duration(0,1,0);
dt_v = startTime:deltaTime:endTime;

%Vents
% centerAlt = 3292;
% center = Location(37+45/60+17.19/3600,14+59/60+48.41/3600,centerAlt); %nec
centerAlt = 3257;
plumeLoc = Location(37+44/60+50.17/3600,15+3.67/3600,3257); %sec
center = copy(plumeLoc);
rcvr = VaprRcvr(37+44/60+50.17/3600,15+3.67/3600,3257,'RCVR');

sim = VaprSim();
[area_lat,area_lon,~] = center.lla;
sim.areaLoc.setLoc(area_lat,area_lon,0);
sim.areaRes = [200,200,200];
sim.setTime(startTime);
sim.addVent(plumeLoc);
sim.minEle = minEle;
sim.addRec(rcvr);

%DEM
sim.loadDEM('dem/Etna_Contour.csv');
% maxNorth = 4000;
% minNorth = -4000;
% maxEast = 4000;
% minEast = -4000;
maxNorth = 5000;
minNorth = -5000;
maxEast = 6000;
minEast = -6000;
east_v = sim.demX(1,:);
north_v = sim.demY(:,1);
east_v = east_v(east_v>minEast & east_v<maxEast);
north_v = north_v(north_v>minNorth & north_v<maxNorth);

%Plume
plumeRad = 200; %m
plumeHeight = 3000; %m
sim.minEle = 0;

%%
%Receivers
sim.minEle = minEle;
vent = plumeLoc;
obsCnt_m = zeros(length(north_v),length(east_v),length(dt_v));
for eIdx = 1:length(east_v)
    t1 = datetime('now');
    eT = east_v(eIdx);
    fprintf('%02.0f:%02.0f:%02.0f',t1.Hour,t1.Minute,t1.Second);
    fprintf(' Receiver Grid: %02.0fx-- Location: %0.0f,--,--\n',...
        eIdx,eT);
    parfor nIdx = 1:length(north_v)
        nT = north_v(nIdx);
        alt(nIdx,eIdx) = sim.findAlt(eT,nT);
        rcvr.setENU(eT,nT,alt(nIdx,eIdx)+antHeight,sim.areaLoc);
        [ obsCntTmp_v, covTemp,~,elSat_v] = ...
            simSiteView( 1, vent, sim, plumeRad, plumeHeight, dt_v );
        obsCnt_m(nIdx,eIdx,:) = obsCntTmp_v;
        obsPos_m(nIdx,eIdx) = covTemp;
        obsEl_m{nIdx}{eIdx} = elSat_v;
        
        el_temp = round(elSat_v);
        el_temp = el_temp(el_temp>0);
        if(~isempty(el_temp))
            obsMeanEl_m(nIdx,eIdx,:) = mean(squeeze(mean(el_temp)));
            obsMinEl_m(nIdx,eIdx,:) = min(squeeze(min(el_temp)));
            obsMaxEl_m(nIdx,eIdx,:) = max(squeeze(max(el_temp)));
        else
        	obsMeanEl_m(nIdx,eIdx,:) = 0;
            obsMinEl_m(nIdx,eIdx,:) = 0;
            obsMaxEl_m(nIdx,eIdx,:) = 0;    
        end
    end
end
[eM,nM] = meshgrid(east_v,north_v);

%%
[mE, mELoc] = max(obsPos_m);
[mNE, mNELoc] = max(mE);
fprintf('Max %0.1f%% %0.0f N %0.0f E\n',mNE,mNELoc,mELoc(mNELoc));

%%
[pENU(1),pENU(2),pENU(3)] = plumeLoc.enu(sim.areaLoc);

epdn = VaprRcvr(37+45/60+57.24/3600,15+1/60+.48/3600,2804,'EPDN');
ecne = VaprRcvr(37+45/60+55.08/3600,15+6.48/3600,2892+50,'ECNE');
eplu = VaprRcvr(37+45/60+54.36/3600,14+59/60+8.52/3600,2913,'EPLU');
ecpn = VaprRcvr(37+44/60+37/3600,14+59/60+11/3600,2991,'ECPN');

rcvr_c = {epdn,ecne,eplu,ecpn};

    
altOut = alt;
%%
figure; hold on;
surf(eM,nM,altOut,obsPos_m,'edgecolor','none');
cl = colorbar;
plot3(pENU(1),pENU(2),pENU(3)+200,'kv',...
    'MarkerSize',10,'markerfacecolor','r','LineWidth',1.5);
% text(pENU(1),pENU(2),pENU(3)+400,' SEC','Color','k',...
%     'HorizontalAlignment','left','FontWeight','bold');
for ii=1:length(rcvr_c)
    [rENU(1),rENU(2),rENU(3)] = rcvr_c{ii}.enu(sim.areaLoc);
    plot3(rENU(1),rENU(2),rENU(3)+50,'wo',...
    'MarkerSize',6,'markerfacecolor','w');
    if(ii<=3)
    text(rENU(1)+100,rENU(2)+100,rENU(3)+500,rcvr_c{ii}.name,'Color','w',...
            'HorizontalAlignment','left','FontSize',8,'FontWeight','Bold',...
            'BackgroundColor','k','Margin',.5);  
    else
    text(rENU(1)-100,rENU(2)-100,rENU(3)+500,rcvr_c{ii}.name,'Color','w',...
            'HorizontalAlignment','right','FontSize',8,'FontWeight','Bold',...
            'BackgroundColor','k','Margin',.5);  
    end
    plot3(rENU(1),rENU(2),rENU(3)+50,'ko',...
    'MarkerSize',6,'markerfacecolor','w','LineWidth',1);
end
caxis([0 100]);
cmap = colormap(jet(90));
cmap(1,:) = [0,0,.4];
cmap(end+1:end+10,:) = ones(10,1)*cmap(end,:);
colormap(cmap);
xlabel('East');
ylabel('North');
zlabel('Altitude');
ylabel(cl,'Coverage (%)');
title(sprintf('Etna %0.1fx%0.1fkm Plume Visibility',plumeRad*2/1000,plumeHeight/1000));
view(30,65);
xlim([minEast,maxEast]);
zlim([500,4000]);
xticks([-5000,0,5000]);
zticks([0,2000,4000]);
mat_fname = sprintf('mat/etna_pos_100_%dx%d.mat',plumeRad,plumeHeight);
% save(mat_fname,'eM','nM','altOut','obsPos_m','obsCnt_m');



%% Mean
figure; hold on;
surf(eM,nM,altOut,obsMeanEl_m,'edgecolor','none');
cl = colorbar;
plot3(pENU(1),pENU(2),pENU(3)+200,'kv',...
    'MarkerSize',10,'markerfacecolor','r','LineWidth',1.5);

for ii=1:length(rcvr_c)
    [rENU(1),rENU(2),rENU(3)] = rcvr_c{ii}.enu(sim.areaLoc);
    plot3(rENU(1),rENU(2),rENU(3)+50,'wo',...
    'MarkerSize',6,'markerfacecolor','w');
    if(ii<=3)
    text(rENU(1)+100,rENU(2)+100,rENU(3)+500,rcvr_c{ii}.name,'Color','w',...
            'HorizontalAlignment','left','FontSize',8,'FontWeight','Bold',...
            'BackgroundColor','k','Margin',.5);  
    else
    text(rENU(1)-100,rENU(2)-100,rENU(3)+500,rcvr_c{ii}.name,'Color','w',...
            'HorizontalAlignment','right','FontSize',8,'FontWeight','Bold',...
            'BackgroundColor','k','Margin',.5);  
    end
    plot3(rENU(1),rENU(2),rENU(3)+50,'ko',...
    'MarkerSize',6,'markerfacecolor','w','LineWidth',1);
end
caxis([20 90]);
cmap = colormap(jet(90));
cmap(1,:) = [0,0,.4];
cmap(end+1:end+10,:) = ones(10,1)*cmap(end,:);
colormap(cmap);
xlabel('East');
ylabel('North');
zlabel('Altitude');
ylabel(cl,'Mean Elevation');
title(sprintf('Etna %0.1fx%0.1fkm Mean Elevation Angle',plumeRad*2/1000,plumeHeight/1000));
view(30,65);
xlim([minEast,maxEast]);
zlim([500,4000]);
xticks([-5000,0,5000]);
zticks([0,2000,4000]);
mat_fname = sprintf('mat/etna_pos_100_%dx%d.mat',plumeRad,plumeHeight);
% save(mat_fname,'eM','nM','altOut','obsPos_m','obsCnt_m');

%% Max
figure; hold on;
surf(eM,nM,altOut,obsMaxEl_m,'edgecolor','none');
cl = colorbar;
plot3(pENU(1),pENU(2),pENU(3)+200,'kv',...
    'MarkerSize',10,'markerfacecolor','r','LineWidth',1.5);

for ii=1:length(rcvr_c)
    [rENU(1),rENU(2),rENU(3)] = rcvr_c{ii}.enu(sim.areaLoc);
    plot3(rENU(1),rENU(2),rENU(3)+50,'wo',...
    'MarkerSize',6,'markerfacecolor','w');
    if(ii<=3)
    text(rENU(1)+100,rENU(2)+100,rENU(3)+500,rcvr_c{ii}.name,'Color','w',...
            'HorizontalAlignment','left','FontSize',8,'FontWeight','Bold',...
            'BackgroundColor','k','Margin',.5);  
    else
    text(rENU(1)-100,rENU(2)-100,rENU(3)+500,rcvr_c{ii}.name,'Color','w',...
            'HorizontalAlignment','right','FontSize',8,'FontWeight','Bold',...
            'BackgroundColor','k','Margin',.5);  
    end
    plot3(rENU(1),rENU(2),rENU(3)+50,'ko',...
    'MarkerSize',6,'markerfacecolor','w','LineWidth',1);
end
caxis([20 90]);
cmap = colormap(jet(90));
cmap(1,:) = [0,0,.4];
cmap(end+1:end+10,:) = ones(10,1)*cmap(end,:);
colormap(cmap);
xlabel('East');
ylabel('North');
zlabel('Altitude');
ylabel(cl,'Maximum Elevation');
title(sprintf('Etna %0.1fx%0.1fkm Maximum Elevation Angle',plumeRad*2/1000,plumeHeight/1000));
view(30,65);
xlim([minEast,maxEast]);
zlim([500,4000]);
xticks([-5000,0,5000]);
zticks([0,2000,4000]);
mat_fname = sprintf('mat/etna_pos_100_%dx%d.mat',plumeRad,plumeHeight);
% save(mat_fname,'eM','nM','altOut','obsPos_m','obsCnt_m');

%%
%EPDN
el_temp = round(obsEl_m{36}{36});
el_temp = el_temp(el_temp>0);
figure;
subplot(3,1,1);
histogram(el_temp);
xlim([19.5,90.5]);
ylim([0 50]);
title('EPDN')
xlabel('Elevation Angle');
ylabel('Minutes');
fprintf('Minutes Visible:%d Mean Elevation:%0.2f\n',sum(el_temp>0),mean(el_temp));

%EPLU
subplot(3,1,2);
el_temp = round(obsEl_m{23}{36});
el_temp = el_temp(el_temp>0);
histogram(el_temp);
xlim([19.5,90.5]);
title('EPLU')
xlabel('Elevation Angle');
ylabel('Minutes');
ylim([0 50]);
fprintf('Minutes Visible:%d Mean Elevation:%0.2f\n',sum(el_temp>0),mean(el_temp));

%ECPN
subplot(3,1,3);
el_temp = round(obsEl_m{23}{24});
el_temp = el_temp(el_temp>0);
histogram(el_temp);
xlim([19.5,90.5]);
title('ECPN')
xlabel('Elevation Angle');
ylabel('Minutes');
ylim([0 50]);
fprintf('Minutes Visible:%d Mean Elevation:%0.2f\n',sum(el_temp>0),mean(el_temp));

%%
clear obsPrc_m;
for ii=1:length(obsCnt_m(:,1,1))
    for jj=1:length(obsCnt_m(1,:,1))
        clear tmp;
        tmp = squeeze(obsCnt_m(ii,jj,:));
        obsPrc_m(ii,jj) = 100*sum(tmp>1)./length(tmp);
    end
end


figure; hold on;
surf(eM,nM,altOut,obsPrc_m,'edgecolor','none');
% surf(eM,nM,altOut);
cl = colorbar;
plot3(pENU(1),pENU(2),pENU(3)+200,'kv',...
    'MarkerSize',10,'markerfacecolor','r','LineWidth',1.5);
for ii=1:length(rcvr_c)
    [rENU(1),rENU(2),rENU(3)] = rcvr_c{ii}.enu(sim.areaLoc);
    plot3(rENU(1),rENU(2),rENU(3)+50,'o',...
    'MarkerSize',6,'markerfacecolor','w');
    text(rENU(1)+100,rENU(2)+100,rENU(3)+500,rcvr_c{ii}.name,'Color','w',...
            'HorizontalAlignment','left','FontSize',8,'FontWeight','Bold',...
            'BackgroundColor','k','Margin',.5);    
    plot3(rENU(1),rENU(2),rENU(3)+50,'ko',...
    'MarkerSize',6,'markerfacecolor','w','LineWidth',1);
end
caxis([0 100]);
cmap = colormap(jet(90));
cmap(1,:) = [0,0,.4];
cmap(end+1:end+10,:) = ones(10,1)*cmap(end,:);
colormap(cmap);
xlabel('East (km)');
ylabel('North (km)');
zlabel('Altitude (km)');
ylabel(cl,'% Time');
title(sprintf('Redoubt %0.1fx%0.1fkm Multiple Plume Intersections',plumeRad*2/1000,plumeHeight/1000));
view(30,75);
zlim([0,4000]);
xlim([minEast,maxEast]);
ylim([minNorth,maxNorth]);
xticks([-5000,0,5000,10000]);
yticks([-5000,0,5000,10000]);
zticks([0,2000,4000]);
xticklabels({'-5','0','5','10'});
yticklabels({'-5','0','5','10'});
zticklabels({'0','2','4'});
