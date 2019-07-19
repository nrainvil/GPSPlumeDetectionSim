clear all;
% close all;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('lib');
addpath('emag_lib');
addpath('plotGen');

minEle = 25; %Degrees
antHeight = 2; % m 

%Simulation time 
startTime = datetime(2013,11,23,12,0,0,'TimeZone','UTC');
endTime = datetime(2013,11,24,11,55,0,'TimeZone','UTC');
deltaTime = duration(0,1,0);
dt_v = startTime:deltaTime:endTime;

%Vents
% centerAlt = 3292;
% center = Location(37+45/60+17.19/3600,14+59/60+48.41/3600,centerAlt); %nec
centerAlt = 3257;
plumeLoc = Location(37+44/60+50.17/3600,15+3.67/3600,3257); %sec
center = copy(plumeLoc);
rcvr = VaprRcvr(37+44/60+50.17/3600,15+3.67/3600,3257,'RCVR');;

sim = VaprSim();
[area_lat,area_lon,~] = center.lla;
sim.areaLoc.setLoc(area_lat,area_lon,0);
sim.setTime(startTime);
sim.addVent(plumeLoc);
sim.minEle = minEle;
sim.addRec(rcvr);

%DEM
sim.loadDEM('dem/Etna_Contour.csv');
maxNorth = 4000;
minNorth = -4000;
maxEast = 4000;
minEast = -4000;
east_v = sim.demX(1,:);
north_v = sim.demY(:,1);
east_v = east_v(east_v>minEast & east_v<maxEast);
north_v = north_v(north_v>minNorth & north_v<maxNorth);

%Plume
plumeRad = 500; %m
plumeHeight = 5000; %m
sim.minEle = 0;
[ ~, ~,az,el] = ...
    simSiteView( 1, plumeLoc, sim, plumeRad, plumeHeight, dt_v );
time_v = dt_v;
prn_v = [1:3,5:32];
% prn_v = [12,19,26];

fig1=figure;
polaraxes(); hold on;
for pIdx=1:length(prn_v)
    prn = prn_v(pIdx);
    azP = (pi()/180).*(az(:,prn));
    elP = (el(:,prn));
    plt{prn} = polarplot(azP,elP,'.');
    set(plt{prn},'MarkerSize',4);
end
for pIdx=1:length(prn_v)
    prn = prn_v(pIdx);
    cl = get(plt{prn},'Color');
    pltT(pIdx) = polarplot(-5,360,'-');
    set(pltT(pIdx),'Color',cl);
    set(pltT(pIdx),'LineWidth',2);
    pTxt{pIdx} = sprintf('%2.0f',prn); 
end

xmax = 120;
ymax = 90;
ysize = 600;
xsize = floor(ysize*(xmax+90)/(ymax+90));
% axis([-90,xmax,-90,ymax]);
pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
pax.RDir = 'reverse';
thetaticks([0,90,180,270]);
thetaticklabels({'N','E','S','W'});
rlim([0 90]);

% yticks([-90, -45, 0, 45, 90]);
lg1 = legend(pltT,pTxt,'Location','northwest');
title(lg1,'PRN');
set(gcf,'Position',[50 50 xsize ysize]);
title('GPS Satellite Position');

% tIdx = find(el(:,12)>0);
% prn = 12;
% hold on;
% for ii=1:10
%     tm = tIdx(1) + (ii-1)*60;
%     polarplot(pi()/180*az(tm,prn),el(tm,prn),'o','Color','b');
% end

az_bin = 45/10;
el_bin = 5;
az_len = floor(360/az_bin);
el_len = floor(90/el_bin);
azelData = zeros(az_len,el_len);
for ii = 1:az_len
    for jj = 1:el_len
        azMax = ii*az_bin;
        elMax = jj*el_bin;
        
        azSl = (az>(azMax-az_bin) & az<=(azMax));
        elSl = (el>(elMax-el_bin) & el<=(elMax));
        azelData(ii,jj) = sum(sum(azSl & elSl));
%         azelData(ii,jj) = sum(sum(azSl & elSl) > 0);
    end
end
% azelData(azelData<15)=0;
fig2 = plot_azel(azelData,az_bin,el_bin,...
                    'GPS Satellite Visibility',...
                    [0 45],'Time Visible (Min)');


set(gcf,'Position',[50 50 xsize ysize]);