% close all;
clear all;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('kml_lib');
addpath('emag_lib');
addpath('plotGen');

dtime=datetime(2013,11,23,10,15,0,'TimeZone','UTC');
etnaSetupAll; %Initialize sim, vents and DEM

%Plume
% Based on
plumeLoc = copy(sec);                      %Bonaccorso
rhoAsh = 2.813*(100^3);                    %g/m^3 Suwanosejima, Oguchi 2009 
% rhoAsh = 1*(100^3);                        %g/m^3 Marzano 2006
% e_s_rel = 5 - 1j*.3;                       %(relative) Adams and Onizuka
e_s_rel = 6 - 1j*.27;                      % Adams
% e_s_rel = 6.109 - 1j*0.136;
% e_s_rel = 85.7-1j*14.1; %Water
% rhoAsh = 9.998e5; %Water

[~,~,plumeBase] = plumeLoc.enu(sim.areaLoc);
plumeHeight = 3000 + plumeBase;            %(m) Lava Fountain
uN = .0164;                                 %(m) Mode of particle size
k = 2.3;
ashDia_v = .0001:.0025:.1001;              % Ash Diamter mm

%Create plume  shifted 150m to the NE
plumeRes = [50,50,40];
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
MER = 3.2e6;                               %kg/s Mass Eruption Rate
Vvent = 230;                               %m/s Vent Velocity
plume.addMEREllCylPlume(sA,sB,dirA,plumeHeight,MER,Vvent,sim.areaLoc);

%
sim.findSignalPaths(10);

for rcvrIdx=1:length(sim.receivers)
    for pathIdx=1:length(sim.signalPaths{rcvrIdx}(1,:))
        ray = sim.signalPaths{rcvrIdx}(pathIdx);
        prn = sim.signalPRNs{rcvrIdx}(pathIdx);
        [pathMass, pathLen, atn] = sim.findSignalAtten(ray);
        if(pathMass>0)
            fprintf('Intersection at %s (%d) ',...
                sim.receivers(rcvrIdx).name,rcvrIdx);
            fprintf('PRN: %d (%d)\n Integrated Mass: %0.0f (kg/m^2),',...
                prn,pathIdx,pathMass*1e-3);
            fprintf(' Length %0.2f (m), Attenuation: %0.4f (dB)\n',...
                pathLen,atn);
            [pathMass, pathLen, atn] = sim.findSignalAtten(ray,1);
        end
    end
end
%%
sim.plotSim();
xlim([-4000,4000]);
ylim([-4000,4000]);
