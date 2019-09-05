warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('kml_lib');
addpath('emag_lib');
clear all;

frame_v = [];
startTime = datetime(2013,11,23,9,50,0,'TimeZone','UTC');
% endTime = datetime(2013,11,23,9,51,0,'TimeZone','UTC');
deltaTime = duration(0,1,0); 
endTime = datetime(2013,11,23,10,17,0,'TimeZone','UTC');

dt_v = startTime:deltaTime:endTime;

for dtIdx = 1:length(dt_v)
    %update plume
    elapsedTime = dt_v(dtIdx)-startTime;
    dt = dt_v(dtIdx);

cfg_m = [6.33-1j*0.204, 2.3, .0164, 3.2e6, 230]; %Nominal

%Receivers
epdn = VaprRcvr(37+45/60+57.24/3600,15+1/60+.48/3600,2804,'EPDN');

%Vents
nec = Location(37+45/60+17.19/3600,14+59/60+48.41/3600,3292);
vor = Location(37+45/60+2.33/3600,14+59/60+36.82/3600,3267);
sec = Location(37+44/60+50.17/3600,15+3.67/3600,3257);

for cIdx=1:length(cfg_m(:,1))
clear vsr sim atn_v
e_s_rel = cfg_m(cIdx,1);
k = cfg_m(cIdx,2);
uN = cfg_m(cIdx,3);
MER = cfg_m(cIdx,4);
Vvent = cfg_m(cIdx,5);

sim = VaprSim();
[area_lat,area_lon,~] = sec.lla;
sim.areaLoc.setLoc(area_lat,area_lon,0);
sim.addRec(epdn);
sim.addVent(nec);
sim.addVent(vor);
sim.addVent(sec);
sim.loadDEM('dem/Etna_Contour.csv');

%Plume
% Based on
plumeLoc = copy(sec);                      %Bonaccorso
rhoAsh = 2.813*(100^3);                    %g/m^3 Suwanosejima, Oguchi 2009 
[~,~,plumeBase] = plumeLoc.enu(sim.areaLoc);
plumeHeight = 3000 + plumeBase;            %(m) Lava Fountain
ashDia_v = .0001:.0025:.1001;              % Ash Diameter mm

%Create plume  shifted 150m to the NE
plumeRes = [40,40,20];
plumeSize = [1000,1000,plumeHeight-plumeBase];
shiftDist = 200;                           %(m) Estimated due to wind
shiftDir = [8,3];                          %NE Direction of Plume Corradini
shiftVec = shiftDist.*shiftDir./norm(shiftDir);
plumeLoc.shiftENU([shiftVec,0],sim.areaLoc);
plume = Plume(plumeLoc,plumeSize,plumeRes);
plume.k = k;
plume.initParticleTypes(ashDia_v,e_s_rel.*ones(size(ashDia_v)),...
                            rhoAsh.*ones(size(ashDia_v)));
sim.addPlume(plume);

%Create Elliptical Cylinder MER plume
plume.uN = uN;
sA = 200;
sB = 100;
dirA = shiftDir;                           
plume.addMEREllCylPlume(sA,sB,dirA,plumeHeight,MER,Vvent,sim.areaLoc);
                
%Simulation time 
sim.setTime(dt);
fig1 = sim.plotSim();

%Plot Formatting
xtick_v = [-4000,-2000,0,2000,4000];
xticks(xtick_v);
xtick_v = xtick_v./1000;
xticklabels({num2str(xtick_v(1)),...
    num2str(xtick_v(2)),...
    num2str(xtick_v(3)),...
    num2str(xtick_v(4)),...
    num2str(xtick_v(5))});

ytick_v = [-4000,-2000,0,2000,4000];
yticks(ytick_v);
ytick_v = ytick_v./1000; %km
yticklabels({num2str(ytick_v),...
    num2str(ytick_v(2)),...
    num2str(ytick_v(3)),...
    num2str(ytick_v(4)),...
    num2str(ytick_v(5))});

xlim([-1000,3000]);
ylim([-1000,3000]);
colormap('bone');
caxis([1000,3300]);
end
frame_v = [frame_v,getframe(fig1)];
end

%% Display Movie
close all;
mfig = figure('pos',[300,100,550,450]);
movie(mfig,frame_v,60,1);