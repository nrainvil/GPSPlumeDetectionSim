warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('kml_lib');
addpath('emag_lib');
profile on;
g = 9.8; %m/s
%Simulation time 
%Event 19
% startTime = datetime(2009,04,04,13,57,0,'TimeZone','UTC');
% endTime = datetime(2009,04,09,14,00,0,'TimeZone','UTC');
% endTime = datetime(2009,04,04,14,23,0,'TimeZone','UTC');

%Event 5
% startTime = datetime(2009,03,23,12,30,0,'TimeZone','UTC');
startTime = datetime(2009,03,23,12,30,0,'TimeZone','UTC');
endTime = datetime(2009,03,23,12,51,0,'TimeZone','UTC');
% endTime = datetime(2009,03,23,12,51,0,'TimeZone','UTC');

% cfg_m = [6.33-1j*0.204, 2.3, .0164, 3.2e6, 230]; %Nominal


% cfg_m = [...
%          6.33-1j*0.204, 2.3, .0164, 3.2e8, 230];% ... %Nominal
%          7.60-1j*0.320, 2.1, .0164, 4.0e6, 220; ...
%          7.60-1j*0.320, 2.1, .0200, 4.0e6, 230; ...
%          85.7-1j*14.1, 2.1, .0164, 4.0e6, 220; ...
%          85.7-1j*14.1, 2.1, .0200, 4.0e6, 230]; % Best


e_s_dryhail = 3.17 - 1j*.004; 
e_s_spghail = 22.6 - 1j*11.41;
e_s_wethail = 80.9 - 1j*23.865;

%e_s_hail, k Etna, uN Etna, MER Redoubt, vVent 14.9k alt 
cfg_m = [e_s_wethail, 2.3, .0164, 1.53e7, sqrt(2*g*14.9e3)]; %Nominal
%Vents
centerAlt = 2.659078673646549e+03;
vent = Location(60.48524,-152.75315,centerAlt);

rvbmAlt = 1.515933747641333e+03 + 2;
rvbm = VaprRcvr(60.487,-152.844,rvbmAlt,'RVBM');
rbed = VaprRcvr(60.454,-152.745,1.490697754927134e+03 + 2,'RBED');
dumm = VaprRcvr(60.580,-152.664,1.725229710264339e+02 + 2,'DUMM');
ac17 = VaprRcvr(60.664,-152.404,6.654326198126273e+02 + 2,'AC17');
% rcvr_c = {rvbm,rbed,dumm};
rcvr_c = {rvbm};%,rbed};
plumeLoc = copy(vent);
rcvr = copy(rvbm);
center = copy(plumeLoc);
           
for cIdx=1:length(cfg_m(:,1))
clear vsr sim atn_v dt_v
e_s_rel = cfg_m(cIdx,1);
k = cfg_m(cIdx,2);
uN = cfg_m(cIdx,3);
MER = cfg_m(cIdx,4);
Vvent = cfg_m(cIdx,5);

sim = VaprSim();
[area_lat,area_lon,~] = center.lla;
sim.areaLoc.setLoc(area_lat,area_lon,0);
sim.setTime(startTime);
sim.addVent(plumeLoc);
sim.minEle = 15;
sim.addRec(rcvr);
sim.addVent(vent);

%DEM
sim.areaSize=[15000,15000,18000];
sim.areaRes = [500,500,500];
sim.loadDEM('dem\redoubt_contour_v4.csv');

%Plume
% Based on
rhoAsh = 2.813*(100^3);                    %g/m^3 Suwanosejima, Oguchi 2009 
[~,~,plumeBase] = plumeLoc.enu(sim.areaLoc);
plumeHeight = 15000 + plumeBase;            %(m) Lava Fountain
ashDia_v = [.001,.002,.004,.008,.016,.031,.063,.125,.250,...
            .500,1,2,4,8,16,32];
ashDia_v = ashDia_v*1e-3; %Ash Diameter m

ashSolidRW = [0,0,0,0,0,0,0,4.14,08.95,...
              11.42,10.27,2.55,3.06,2.67,0.83,0.04]; % Weight% single part.
ashAggRW =   [0,0,0,0,0,0,0,5.30,21.32,...
              24.63,02.86,1.96,0.01,0.00,0.00,0.00]; % Weight% Aggregates
          
%Convert % to kg/m^3
% typeRho_v = 

%Create plume Time, sA, sB, shiftDist, shiftDir E,N, height
%Constant density - same MER results in different densities
cfg_time_m =    ...
    [datetime(2009,03,23,12,30,0,'TimeZone','UTC');...
     datetime(2009,03,23,12,33,0,'TimeZone','UTC'); ...
     datetime(2009,03,23,12,41,0,'TimeZone','UTC'); ...
     datetime(2009,03,23,12,50,0,'TimeZone','UTC')];
cfg_plume_m = ...
    [0,0,0,-.1,1,plumeHeight,MER,Vvent;...
     3000,2500,2500,-.1,1,plumeHeight,MER,Vvent; ...
     6500,6000,5500,-.1,1,plumeHeight,MER,Vvent; ...
     100,100,5500,-.1,1,plumeHeight,MER,Vvent];

plumeRes = [140,140,140];
%plumeRes = [500,500,500];
plumeSize = [14000,14000,plumeHeight-plumeBase];
plume = Plume(plumeLoc,plumeSize,plumeRes);
plume.k = k;
plume.initParticleTypes(ashDia_v,e_s_rel.*ones(size(ashDia_v)),...
                            rhoAsh.*ones(size(ashDia_v)));
sim.addPlume(plume);

%Create Elliptical Cylinder MER plume
plume.uN = uN;                    

%Run sim and create output plots
gifFile = sprintf('plots/redoubt_evolv_140.gif');

vsr = VaprSimRun(startTime,endTime,sim);
vsr.deltaTime = duration(0,1,0);
vsr.cfg_time = cfg_time_m;
vsr.cfg_plume = cfg_plume_m;

%Run sim and create output plots
createGif = 1;
%vsr.sim.cfg_plume = cfg_plume_m;
vsr.sim.plotLim = [-12500,12500];
vsr.sim.plotView = [0,90]; 
t1 = datetime('now');
vsr.runSim(createGif); %Run Sim

t2 = datetime('now');
fprintf('Start Time: %02.0f:%02.0f:%02.0f \n',t1.Hour,t1.Minute,t1.Second);
fprintf('End Time: %02.0f:%02.0f:%02.0f \n',t2.Hour,t2.Minute,t2.Second);

vsr.plotSim();
matFile = sprintf('mat/rdbt_l2_%2.0f_%02.0f_%02.0f_%03.0f_%03.0f_%03.0f.mat',...
                            real(e_s_rel)*100,abs(imag(e_s_rel))*100,...
                            k*10,uN*1000,MER*1e-5,Vvent);

dt_v = vsr.startTime:vsr.deltaTime:vsr.sim.time;
for dt=1:length(vsr.atn_c)
    atn_v(dt) = max(abs(vsr.atn_c{dt}{1}));
end
% save(matFile,'dt_v','atn_v');
end

if(createGif==1)
    mfig = figure('pos',[300,100,950,950]);
    fprintf('Saving %s\n',gifFile);
    if exist(gifFile, 'file')==2
        delete(gifFile);
    end
    for ii=1:length(vsr.frame_v)
        im = frame2im(vsr.frame_v(ii));
        [imIdx,cm] = rgb2ind(im,256);
        if(ii==1)
            imwrite(imIdx,cm,gifFile,'gif','DelayTime',1,'Loopcount',inf);
        else
            imwrite(imIdx,cm,gifFile,'gif','DelayTime',1,'WriteMode','append');
        end
    end
    movie(mfig,vsr.frame_v,60,1);
end