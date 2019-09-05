clear all;
% close all;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
addpath('kml_lib');
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
rcvr = VaprRcvr(37+44/60+50.17/3600,15+3.67/3600,3257,'RCVR');

sim = VaprSim();
[area_lat,area_lon,~] = center.lla;
sim.areaLoc.setLoc(area_lat,area_lon,0);
sim.areaRes = [250,250,250];
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
plumeWidth = 2*plumeRad;
sim.minEle = 0;

%Receivers
[eM,nM] = meshgrid(east_v,north_v);
mat_fname = sprintf('mat/etna_pos_%dx%d.mat',plumeRad,plumeHeight);
load(mat_fname);
%%
obsCntOrig_m = obsCnt_m;
obsPosOrig_m = obsPos_m;
sim.minEle = minEle;
vent = plumeLoc;
eMax = length(east_v)-1;
nMax = length(north_v)-1;

radLen = 2;
tanLen = 8;
rMin_v = [1000,2000,3000];
rMax_v = [1300,2300,3300];

obsCntTotAll = zeros(length(obsCnt_m(1,1,:)),1);
for jj=1:radLen
rMin = rMin_v(jj);
rMax = rMax_v(jj);
obsCnt_m = obsCntOrig_m;
obsPos_m = obsPosOrig_m;
for eIdx = 1:length(east_v)
    eT = east_v(eIdx);
    for nIdx = 1:length(north_v)
        nT = north_v(nIdx);
        alt = altOut(nIdx,eIdx);
        if(norm([eT,nT])<rMin || norm([eT,nT])>rMax || ...
           eIdx>eMax || nIdx>nMax)
            obsPos_m(nIdx,eIdx)=0;
            obsCnt_m(nIdx,eIdx,:) = zeros(length(dt_v),1);
        end
    end
end

[pENU(1),pENU(2),pENU(3)] = plumeLoc.enu(sim.areaLoc);
obsCntTot = zeros(length(obsCnt_m(1,1,:)),1);
oCNew = zeros(length(north_v),length(east_v),length(dt_v));
cov = 0;
covMin = 0;
if(jj==1)
    tanLenNew = tanLen-1;
else
    tanLenNew = tanLen;
end
for ii=1:tanLenNew
    %Max Location and fills gap
    for eIdx = 1:eMax
        for nIdx = 1:nMax
            oCNew(nIdx,eIdx,:) = obsCntTot;
        end
    end
    
    oCOr = or(obsCnt_m>0,oCNew>covMin);
    if(cov<1)
        [mE,mELoc] = max(100*squeeze(sum(oCOr,3)./length(dt_v))); %??
    else
        [mE, mELoc] = max(obsPos_m);
    end
  
    [mNE, mNELoc] = max(mE);
    fprintf('%0.0f: %0.1f%% %0.0f N %0.0f E\n',ii,mNE,mELoc(mNELoc),mNELoc);
    nIdx1 = mELoc(mNELoc);
    eIdx1 = mNELoc;

    eT1 = east_v(eIdx1)+sim.areaRes(1)/2;
    nT1 = north_v(nIdx1)+sim.areaRes(2)/2;
    r1ENU = [eT1,nT1,altOut(nIdx1,eIdx1)];
    rnENU(jj,ii,:) = r1ENU;
    obsCntTot = obsCntTot + squeeze(obsCnt_m(nIdx1,eIdx1,:));
    
    for eIdx = 1:eMax
    eT = east_v(eIdx)+sim.areaRes(1)/2;
    for nIdx = 1:nMax
        nT = north_v(nIdx)+sim.areaRes(2)/2;
        r2ENU = [eT,nT,altOut(nIdx,eIdx)];
        %Find tangential distance 
        Rpr1 = r1ENU(1:2) - pENU(1:2);
        Rr1r2 = r2ENU(1:2) - r1ENU(1:2);
        norm(Rr1r2);
        RProj = dot(Rr1r2,Rpr1/norm(Rpr1));
        Rt = sqrt(norm(Rr1r2)^2 - norm(RProj)^2);
        THr = acosd(dot(Rpr1./norm(Rpr1),Rr1r2./norm(Rr1r2)));
        if(Rt<plumeWidth && (norm(Rpr1)>norm(Rr1r2) || THr < 90))
            obsPos_m(nIdx,eIdx)=0;
            obsCnt_m(nIdx,eIdx,:) = zeros(length(dt_v),1);
        end
    end
    end
    if(cov<1)
        cov = sum(obsCntTot>covMin)./length(dt_v);
    end
end
    obsCntTotAll = obsCntTotAll + obsCntTot;
end


%%
% cov = sum(obsCntTot>0)./length(dt_v);
cov = sum(obsCntTotAll>0)./length(dt_v);
figure;
bar(dt_v,obsCntTotAll,'BarWidth',1);
xlim([min(dt_v),max(dt_v)]);
fprintf('%0.0f Sites: %0.2f%% Coverage\n',tanLen*radLen,100*cov);
ylim([0 15]);

%%
stnClr_c = {'g',[1,0,1],[1,.5,0]};

figure; hold on;
surf(eM,nM,altOut,obsPosOrig_m);
% surf(eM,nM,altOut);
cl = colorbar;
plot3(pENU(1),pENU(2),pENU(3)+200,'kv',...
    'MarkerSize',10,'markerfacecolor','r','LineWidth',1.5);
% for ii=1:length(rcvr_c)
%     [rENU(1),rENU(2),rENU(3)] = rcvr_c{ii}.enu(sim.areaLoc);
%     text(rENU(1)+100,rENU(2)+100,rENU(3)+500,rcvr_c{ii}.name,'Color','w',...
%             'HorizontalAlignment','left','FontSize',8,'FontWeight','Bold',...
%             'BackgroundColor','k','Margin',.5);    
%     plot3(rENU(1),rENU(2),rENU(3)+50,'ko',...
%     'MarkerSize',6,'markerfacecolor','w','LineWidth',1);
% end

%Plot new network
for jj=1:radLen
    stnClr = stnClr_c{jj};
for ii=1:tanLen
    plot3(rnENU(jj,ii,1),rnENU(jj,ii,2),rnENU(jj,ii,3)+200,'ko',...
    'MarkerSize',6,'markerfacecolor',stnClr,'LineWidth',1);
    text(rnENU(jj,ii,1)+150,rnENU(jj,ii,2)+150,rnENU(jj,ii,3)+500,...
            sprintf('%0.0f',(jj-1)*tanLen+ii),'Color','w',...
            'HorizontalAlignment','left','FontSize',8,'FontWeight','Bold',...
            'BackgroundColor','k','Margin',.5);    
end
end

caxis([0 100]);
cmap = colormap(jet(90));
cmap(1,:) = [0,0,.4];
cmap(end+1:end+10,:) = ones(10,1)*cmap(end,:);
colormap(cmap);
xlabel('East (km)');
ylabel('North (km)');
zlabel('Altitude (km)');
ylabel(cl,'Coverage (%)');
title(sprintf('Etna %0.1fx%0.1fkm Plume Visibility',plumeRad*2/1000,plumeHeight/1000));
% view(30,65);
view(0,90);
zlim([0,4000]);
xlim([minEast,maxEast]);
ylim([minNorth,maxNorth]);
xticks([-5000,0,5000,10000]);
yticks([-5000,0,5000,10000]);
zticks([0,2000,4000]);
xticklabels({'-5','0','5','10'});
yticklabels({'-5','0','5','10'});
zticklabels({'0','2','4'});

mat_fname = sprintf('mat/etna_net_%dx%d.mat',radLen,tanLen);
save(mat_fname,'dt_v','obsCntTotAll','cov');
%%
clear all;
startTime = datetime(2013,11,23,12,0,0,'TimeZone','UTC');
endTime = datetime(2013,11,24,11,55,0,'TimeZone','UTC');
deltaTime = duration(0,1,0);
dt_v = startTime:deltaTime:endTime;
tanLen = 8;

figure;
hold on;
stnClr_c = {'g',[1,0,1],[1,.5,0]};
for ii = 1:2
    radLen = 3-ii;
    stnClr = stnClr_c{radLen};
    mat_fname = sprintf('mat/etna_net_%dx%d.mat',radLen,tanLen);
    net_c{radLen} = load(mat_fname);


    br_m(radLen) = bar(dt_v,net_c{radLen}.obsCntTotAll,'BarWidth',1,'FaceColor',stnClr);
    xlim([min(dt_v),max(dt_v)]);
    ylim([0 15]);
    xlabel('Time');
    ylabel('Signal Crossings')
end
legend(br_m,'1km','1km,2km');
title('Etna Network Plume Coverage');
