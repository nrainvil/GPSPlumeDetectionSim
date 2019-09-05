clear all;
% close all;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('kml_lib');
addpath('emag_lib');
addpath('plotGen');

dtime=datetime(2013,11,23,0,0,0,'TimeZone','UTC');
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
sB = 200;
dirA = shiftDir;                           
MER = 3.2e6;                               %kg/s Mass Eruption Rate
Vvent = 230;                               %m/s Vent Velocity
plume.addMEREllCylPlume(sA,sB,dirA,plumeHeight,MER,Vvent,sim.areaLoc);

                
%Simulation time 
startTime = datetime(2013,11,23,14,22,0,'TimeZone','UTC');
endTime = datetime(2013,11,23,14,22,0,'TimeZone','UTC');

%Run sim and create output plots
sim.setTime(startTime);
sim.findSignalPaths(20);
sim.plotSim();

createGif = 0;
vsr = VaprSimRun(startTime,endTime,sim);
vsr.deltaTime = duration(0,1,0);
vsr.runSim(createGif);
% vsr.plotSim();

dt_v = vsr.startTime:vsr.deltaTime:vsr.sim.time;
rcvrLen = length(vsr.sim.receivers);
dtLen = length(dt_v);
atn_m = zeros(dtLen,rcvrLen,32);
txt_c = {};

for dtIdx=1:length(dt_v)
    txt_c{dtIdx} = {};
    for rcvrIdx=1:length(vsr.sim.receivers)
        rcvrName = vsr.sim.receivers(rcvrIdx).name;
        pathLen = length(vsr.paths_c{dtIdx}{rcvrIdx});
        for pathIdx=1:pathLen
            atn_m(dtIdx,rcvrIdx,pathIdx) = vsr.atn_c{dtIdx}{rcvrIdx}(pathIdx);
            atn = atn_m(dtIdx,rcvrIdx,:);
            atn(atn>0) = 1;
            if(atn(pathIdx)>0)
                prn = vsr.prns_c{dtIdx}{rcvrIdx}(pathIdx);
                txt_c{dtIdx}{end+1} = {rcvrName,prn};
            end
        end
    end
end

%%
cov=0;
obsCnt_v = zeros(size(dt_v));
for dtIdx=1:length(dt_v)
   dt = dt_v(dtIdx);
   if(~isempty(txt_c{dtIdx})>0)
        fprintf('%02.0f:%02.0f',dt.Hour,dt.Minute);
        txtLen = length(txt_c{dtIdx});
        obsCnt_v(dtIdx) = txtLen;
        for txtIdx=1:txtLen;
            fprintf(' %s:%d',txt_c{dtIdx}{txtIdx}{1},txt_c{dtIdx}{txtIdx}{2});
        end
        cov = cov + 1;
        fprintf('\n');
   end
end

covPerc = 100*cov/length(dt_v);
fprintf('%0.2f%% Coverage\n',covPerc);

%%
figure;
bar(dt_v,obsCnt_v);
yticks([0,1,2,3,4]);
