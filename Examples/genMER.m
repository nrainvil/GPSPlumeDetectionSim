warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('kml_lib');
addpath('emag_lib');
clear all;
 
%Simulation time 
startTime = datetime(2013,11,23,9,59,0,'TimeZone','UTC');
endTime = datetime(2013,11,23,10,00,0,'TimeZone','UTC');
% startTime = datetime(2013,11,23,9,50,0,'TimeZone','UTC')
%endTime = datetime(2013,11,23,10,17,0,'TimeZone','UTC');

%Configuration Matrix
%Row: e_s_rel, k, uN, MER, vVent
cfg_m = [];
cfg_m = [6.33-1j*0.204, 2.3, .0164, 3.2e6, 258]; %Nominal

% cfg_m = [cfg_m;...
%          6.33-1j*0.204, 2.3, .0129, 3.2e6, 230; ...% u_n
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 230; ...
%          6.33-1j*0.204, 2.3, .0200, 3.2e6, 230; ...
%          6.33-1j*0.204, 2.3, .0225, 3.2e6, 230; ...
%          6.33-1j*0.204, 2.3, .0250, 3.2e6, 230];
%      
% cfg_m = [cfg_m;...
%          6.33-1j*0.204, 2.3, .0164, 2.7e6, 230; ...% MER
%          6.33-1j*0.204, 2.3, .0164, 3.0e6, 230; ...    
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 230; ...     
%          6.33-1j*0.204, 2.3, .0164, 3.5e6, 230; ...   
%          6.33-1j*0.204, 2.3, .0164, 3.75e6, 230; ...   
%          6.33-1j*0.204, 2.3, .0164, 4.0e6, 230];    
% 
% cfg_m = [cfg_m;...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 150; ...% vVent   
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 155; ... 
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 160; ...   
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 165; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 170; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 175; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 180; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 185; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 190; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 195; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 200; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 205; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 210; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 215; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 220; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 225; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 250; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 275; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 300; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 325; ...

% cfg_m = [cfg_m;...         
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 350; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 378]; %
% 
% cfg_m = [cfg_m;...
%          7.60-1j*0.320, 2.1, .0164, 4.0e6, 220; ...
%          7.60-1j*0.320, 2.1, .0200, 4.0e6, 230; ...
%          85.7-1j*14.1, 2.1, .0164, 4.0e6, 220; ...
%          85.7-1j*14.1, 2.1, .0200, 4.0e6, 230]; % Best
     
     
% cfg_m = [cfg_m;...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 230; ...%Nominal
%          6.33-1j*0.090, 2.3, .0164, 3.2e6, 230; ...%Min eps''
%          6.33-1j*0.320, 2.3, .0164, 3.2e6, 230; ...%Max eps''
%          5.06-1j*0.090, 2.3, .0164, 3.2e6, 230; ...%Min eps' Min eps''
%          5.06-1j*0.204, 2.3, .0164, 3.2e6, 230; ...%Min eps'
%          5.06-1j*0.320, 2.3, .0164, 3.2e6, 230; ...%Min eps' Max eps''
%          7.60-1j*0.090, 2.3, .0164, 3.2e6, 230; ...%Min eps' Min eps''
%          7.60-1j*0.204, 2.3, .0164, 3.2e6, 230; ...%Max eps'
%          7.60-1j*0.320, 2.3, .0164, 3.2e6, 230; ...%Min eps' Max eps''
%          85.7-1j*14.1, 2.3, .0164, 3.2e6, 230; ...%H20
%          31.5-1j*196, 2.3, .0164, 3.2e6, 230];   %H2SO4
     
% cfg_m = [cfg_m;...
%          6.33-1j*0.204, 2.0, .0164, 3.2e6, 230; ...% k
%          6.33-1j*0.204, 2.1, .0164, 3.2e6, 230; ...
%          6.33-1j*0.204, 2.2, .0164, 3.2e6, 230; ...
%          6.33-1j*0.204, 2.3, .0164, 3.2e6, 230; ...
%          6.33-1j*0.204, 2.4, .0164, 3.2e6, 230; ...
%          6.33-1j*0.204, 2.5, .0164, 3.2e6, 230];

% for cIdx=1:length(cfg_m(:,1))
%     cfg_new(cIdx,:) = [cfg_m(cIdx,:),0];
%     cfg_new(4*(cIdx-1)+1,:) = [cfg_m(cIdx,:),0];
%     cfg_new(4*(cIdx-1)+2,:) = [cfg_m(cIdx,:),1e-6];
%     cfg_new(4*(cIdx-1)+3,:) = [cfg_m(cIdx,:),1e-5];
%     cfg_new(4*(cIdx-1)+4,:) = [cfg_m(cIdx,:),2.7e-5]; 
% end
% cfg_m = cfg_new;
     
% global g_cm_csv
% if(isempty(g_cm_csv))
%     g_cm_csv = csvread('vsr_chrg.csv');
% end

%Longer cfg_m for surface charge
% cfg_m = [...
%         6.33-1j*0.204, 2.3, .0164, 3.2e6, 230,1e-8; ...
%         6.33-1j*0.204, 2.3, .0164, 3.2e6, 230,1e-7; ...
%         6.33-1j*0.204, 2.3, .0164, 3.2e6, 230,1e-6; ...
%         6.33-1j*0.204, 2.3, .0164, 3.2e6, 230,1e-5; ...
%         6.33-1j*0.204, 2.3, .0164, 3.2e6, 230,2.7e-5];

%%     

%Receivers
epdn = VaprRcvr(37+45/60+57.24/3600,15+1/60+.48/3600,2804,'EPDN');

%Vents
nec = Location(37+45/60+17.19/3600,14+59/60+48.41/3600,3292);
vor = Location(37+45/60+2.33/3600,14+59/60+36.82/3600,3267);
sec = Location(37+44/60+50.17/3600,15+3.67/3600,3257);

%Iterate over each entry in cfg_m
for cIdx=1:length(cfg_m(:,1))
clear vsr sim atn_v dt_v

sim = VaprSim();

e_s_rel = cfg_m(cIdx,1);
k = cfg_m(cIdx,2);
uN = cfg_m(cIdx,3);
MER = cfg_m(cIdx,4);
Vvent = cfg_m(cIdx,5);
%For charged surface sim
if(length(cfg_m(1,:))>5)
    eta = cfg_m(cIdx,6);
    sim.eta = eta;
    sim.mie = 2;
end

[area_lat,area_lon,~] = sec.lla;
sim.areaLoc.setLoc(area_lat,area_lon,0);
sim.addRec(epdn);
sim.addVent(nec);
sim.addVent(vor);
sim.addVent(sec);
sim.loadDEM('dem/Etna_Contour.csv');

l1_freq = 1.57542e9; %L1 GPS frequency (Hz)
l2_freq = 1.22760e9; %L1 GPS frequency (Hz)
sim.freq = l1_freq;                  %Hz L1 !!!!!

%Plume
% Based on
plumeLoc = copy(sec);                      %Bonaccorso
rhoAsh = 2.813*(100^3);                    %g/m^3 Suwanosejima, Oguchi 2009 
[~,~,plumeBase] = plumeLoc.enu(sim.areaLoc);
plumeHeight = 8000 + plumeBase;            %(m) Lava Fountain
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

%Run sim and create output plots
gifFile = sprintf('plots/epdn_112313_un%03d_%d%ddir_%dx%d.gif',...
                    1e3*uN,dirA(1),dirA(2),sA,sB);

vsr = VaprSimRun(startTime,endTime,sim);
vsr.deltaTime = duration(0,1,0);

%Run sim and create output plots
createGif = 1;
t1 = datetime('now');
vsr.runSim(createGif);
t2 = datetime('now');
fprintf('Start Time: %02.0f:%02.0f:%02.0f \n',t1.Hour,t1.Minute,t1.Second);
fprintf('End Time: %02.0f:%02.0f:%02.0f \n',t2.Hour,t2.Minute,t2.Second);
vsr.plotSim();

if(createGif==1)
   mfig = figure('pos',[300,100,550,450]);
   movie(mfig,vsr.frame_v,60,1); 
end