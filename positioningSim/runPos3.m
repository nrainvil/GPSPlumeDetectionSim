clear all;
% close all;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('kml_lib');
addpath('emag_lib');
addpath('plotGen');

minEle = 25; %Degrees

%Simulation time 
startTime = datetime(2013,11,23,0,0,0,'TimeZone','UTC');
endTime = datetime(2013,11,23,23,55,0,'TimeZone','UTC');

%Receivers
epdn = VaprRcvr(37+45/60+57.24/3600,15+1/60+.48/3600,2804,'EPDN');
ecne = VaprRcvr(37+45/60+55.08/3600,15+6.48/3600,2892,'ECNE');
eplu = VaprRcvr(37+45/60+54.36/3600,14+59/60+8.52/3600,2913,'EPLU');
ecpn = VaprRcvr(37+44/60+37/3600,14+59/60+11/3600,2991,'ECPN');

%Vents
nec = Location(37+45/60+17.19/3600,14+59/60+48.41/3600,3292);
vor = Location(37+45/60+2.33/3600,14+59/60+36.82/3600,3267);
sec = Location(37+44/60+50.17/3600,15+3.67/3600,3257);

sim = VaprSim();
[area_lat,area_lon,~] = sec.lla;
sim.areaLoc.setLoc(area_lat,area_lon,0);
sim.setTime(startTime);
sim.addRec(epdn);
sim.addRec(ecne);
sim.addRec(eplu);
sim.addRec(ecpn);

% sim.addVent(nec);
% sim.addVent(vor);
sim.addVent(sec)
sim.minEle = minEle;

%DEM
sim.loadDEM('dem/Etna_Contour.csv');

%Plume
% Based on
plumeLoc = copy(sec);                      %Bonaccorso

%Create plume  shifted 150m to the NE
shiftDist = 200;                           %(m) Estimated due to wind
shiftDir = [8,3];                          %NE Direction of Plume Corradini
shiftVec = shiftDist.*shiftDir./norm(shiftDir);
plumeLoc.shiftENU([shiftVec,0],sim.areaLoc);

deltaTime = duration(0,1,0);
dt_v = startTime:deltaTime:endTime;

rcvrLen = length(sim.receivers);
ventLen = length(sim.vents);
dtLen = length(dt_v);
atn_m = zeros(dtLen,rcvrLen,32);
txt_c = {};

plumeRad = 200; %m
plumeHeight = 2716.21; %m

plotOn = 1;
vIdx = 1;
frame_v = [];
obsCnt_v = zeros(size(dt_v));
for rIdx=1:rcvrLen
    fprintf('Receiver: %s\n',sim.receivers(rIdx).name);
    for vIdx=1:ventLen
        vent = plumeLoc;
        [ obsCntTmp_v, covTemp ] = ...
            simSiteView( rIdx, vent, sim, plumeRad, plumeHeight, dt_v );
        obsCnt_v = obsCnt_v + obsCntTmp_v;
    end
end

cov = 0;
for dtIdx=1:length(dt_v)
    if(obsCnt_v(dtIdx)>0)
        cov = cov + 1;
    end
end

covPerc = 100*cov/length(dt_v);
fprintf('%0.2f%% Coverage\n',covPerc);

%%
figure;
bar(dt_v,obsCnt_v);
yticks([0,1,2,3,4]);
xlim([startTime,endTime]);
