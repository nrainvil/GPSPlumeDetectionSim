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
startTime = datetime(2013,11,23,0,0,0,'TimeZone','UTC');
endTime = datetime(2013,11,23,3,55,0,'TimeZone','UTC');

%Run sim and create output plots
createGif = 0;

deltaTime = duration(0,1,0);
dt_v = startTime:deltaTime:endTime;
rcvrLen = length(sim.receivers);
dtLen = length(dt_v);
atn_m = zeros(dtLen,rcvrLen,32);
txt_c = {};

plume_rad = 200; %m
plume_height = 2716.21; %m
center = sim.areaLoc;
ele=25;
plotOn = 1;

vIdx = 1;
frame_v = [];
for dtIdx=1:length(dt_v)
    dt = dt_v(dtIdx);
    sim.setTime(dt);
    sim.findSignalPaths(ele);
    txt_c{dtIdx} = {};
    if(plotOn)
        fig1 = figure('Visible','off','pos',[100,100,900,900]);
        hold on;
    end
    for rcvrIdx=1:length(sim.receivers)
        rcvr = sim.receivers(rcvrIdx);
        [rENU(1), rENU(2), rENU(3)] = rcvr.enu(center);
        for ventIdx=1:length(sim.vents)
            vent = plumeLoc;
            %Plot Vents
            [vENU(1), vENU(2), vENU(3)] = vent.enu(center);
            if(plotOn)
            plot(vENU(1),vENU(2),'dr','LineWidth',2);
            %Plot AZ Limits
            viscircles([vENU(1),vENU(2)],plume_rad,'LineStyle',':');
            end
            for pIdx=1:length(sim.signalPaths{rcvrIdx})
                path = sim.signalPaths{rcvrIdx}(pIdx);
                prn = sim.signalPRNs{rcvrIdx}(pIdx);
                [azel(1),azel(2)] = path.getAzEl();
                azel(1) = gpsAz(azel(1));
%                 fprintf('%s PRN:%d Az: %0.2f El: %0.2f\n',...
%                     rcvr.name,prn,azel(1),azel(2));
                [edge_dist, delta_height, ppENU ] = ...
                    findPP( vent, rcvr, azel, plume_rad, plume_height, center);
                if(edge_dist>0)
                    txt_c{dtIdx}{end+1} = {rcvr.name,prn};
                    az_adj = matAz(azel(1));
                    az_ref_v = edge_dist*[cosd(az_adj),sind(az_adj)]; %Direction vector
                    fprintf('%s intersect on PRN:%d',rcvr.name,prn);
                    fprintf(' Dist: %0.2f m Height: %0.2f m\n',...
                        edge_dist,delta_height);
                    if(plotOn)
                    line([rENU(1),rENU(1)+az_ref_v(1)],...
                        [rENU(2),rENU(2)+az_ref_v(2)],...
                        'LineWidth',2,'LineStyle','-','Color',[0 1 0]);
                    textENU = [rENU(1)+.5*az_ref_v(1),rENU(2)+.5*az_ref_v(2)];
                    text(textENU(1)+75,textENU(2)-75,sprintf('%d',prn),...
                            'Color',[0 1 0]);
                    end
                end
            end
            if(plotOn)
            set(gca,'Color',[.9,.9,.9])
            [vENU(1), vENU(2), vENU(3)] = vent.enu(center);
            [az,vTang] = findAzLim(vENU, rENU, plume_rad);
            vent_az_c{1} = 4000*[cosd(az(1)),sind(az(1)); ...
                cosd(az(2)),sind(az(2))];
            az_c{rcvrIdx} = az;
            end
        end
    end

    if(plotOn)
        for rIdx=1:length(sim.receivers)
            rcvr = sim.receivers(rIdx);
            [rENU(1), rENU(2), rENU(3)] = rcvr.enu(center);
            %Plot Max AZ Limit
            az = az_c{rIdx};
            az_max = max(max(az));
            az_min = min(min(az));
            az_max_v = 10000*[cosd(az_max),sind(az_max)];
            az_min_v = 10000*[cosd(az_min),sind(az_min)];
            line([rENU(1),rENU(1)+az_max_v(1)],[rENU(2),rENU(2)+az_max_v(2)],...
                'LineWidth',.75);
            line([rENU(1),rENU(1)+az_min_v(1)],[rENU(2),rENU(2)+az_min_v(2)],...
                'LineWidth',.75);
            
            %Plot Receivers
            plot(rENU(1),rENU(2),'ok','LineWidth',2,'MarkerSize',10);
            text(rENU(1)+75,rENU(2)+75,rcvr.name);
        end
        xlim([-2000,3000]);
        ylim([-2000,3000]);    
        title(sprintf('%02.0f:%02.0f',dt.Hour,dt.Minute));
        frame_v = [frame_v;getframe(fig1)];
        close(fig1);
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
        for txtIdx=1:txtLen
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

%%
if(plotOn)
mfig = figure('pos',[300,100,900,900]);
movie(mfig,frame_v,1,12);
end
