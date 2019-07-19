function [ fig ] = plotSim( obj, hideFig )
%PLOTSIM 3D Plot of Simulation Environment
%Input: obj (VaprSim)
%Output: fig (Figure)

if(nargin<=1)
    hideFig=0;
end

min_z = 0;

if(hideFig==1)
    fig = figure('Visible','off','pos',[100,100,900,900]);
else
    fig = figure;
end

hold on;
%Plot Receiver Locations and signal paths
for rcvrIdx=1:length(obj.receivers)
    rcvr = obj.receivers(rcvrIdx);
    [rcvrENU(1),rcvrENU(2),rcvrENU(3)] = ...
        rcvr.enu(obj.areaLoc);
    if(rcvrENU(3)<min_z)
        min_z = rcvrENU(3);
    end
    %Receiver Location
    plot3(rcvrENU(1),rcvrENU(2),rcvrENU(3),'ko',...
        'MarkerSize',8,'markerfacecolor','k');
    
    %Signal Paths
%     for pathIdx=1:length(obj.signalPaths{rcvrIdx}(1,:))
%         ray = obj.signalPaths{rcvrIdx}(pathIdx);
%         rayBnd = obj.signalBoundary(ray);
%         x_l = [rcvrENU(1) rayBnd(1)];
%         y_l = [rcvrENU(2) rayBnd(2)];
%         z_l = [rcvrENU(3) rayBnd(3)];
%         line(x_l,y_l,z_l,'Color','g');
%         
%     end
end

%Limits
x_limit = [-1*obj.areaSize(1)/2, obj.areaSize(1)/2];
y_limit = [-1*obj.areaSize(2)/2, obj.areaSize(2)/2];
z_limit = [min_z, min_z+obj.areaSize(3)];

% xlim(x_limit);
% ylim(y_limit);
zlim(z_limit);
xlim(obj.plotLim);
ylim(obj.plotLim);
view(obj.plotView(1),obj.plotView(2));

%Plot DEM
surf(obj.demX,obj.demY,obj.demZ,'FaceColor','interp',...
    'LineStyle','none','FaceLighting','gouraud');

%Plot Vents
for vent_idx=1:length(obj.vents)
    vent = obj.vents(vent_idx);
    [vent_enu(1),vent_enu(2),vent_enu(3)] = ...
        vent.enu(obj.areaLoc);
    plot3(vent_enu(1),vent_enu(2),vent_enu(3),'r^',...
        'MarkerSize',6,'markerfacecolor','r');
end

%Plot Center Location
plot3(0,0,0,'+','Color',[.9,.9,.9],'MarkerSize',4,...
    'markerfacecolor',[.9,.9,.9]);

%Plot Signals
for rcvrIdx=1:length(obj.receivers)
    for pathIdx=1:length(obj.signalPaths{rcvrIdx}(1,:))
        %Signal Paths
        ray = obj.signalPaths{rcvrIdx}(pathIdx);
        rayOrigin = ray.origin.*obj.areaRes;
        rayBnd = obj.signalBoundary(ray);
        x_l = [rayOrigin(1) rayBnd(1)];
        y_l = [rayOrigin(2) rayBnd(2)];
        z_l = [rayOrigin(3) rayBnd(3)];
        
        [pathMass, pathLen, atn] = obj.findSignalAtten(ray);
        if(pathMass>0)
            atnWidth = atn;
            atnWidth(isnan(atnWidth))=1;
            line(x_l,y_l,z_l,'Color','g','LineWidth',atnWidth);
        else
            if(obj.signalPRNs{rcvrIdx}(pathIdx)==4)
              line(x_l,y_l,z_l,'Color','m');
            elseif(obj.signalPRNs{rcvrIdx}(pathIdx)==4)
              line(x_l,y_l,z_l,'Color','c');
            else
              line(x_l,y_l,z_l,'Color','r');
            end
        end
    end
end
    
%Plot Plumes
ver = [1 1 0;0 1 0;0 1 1;...
    1 1 1;0 0 1;1 0 1;
    1 0 0;0 0 0];

fac = [1 2 3 4;4 3 5 6;6 7 8 5;
    1 2 8 7;6 7 1 4;2 3 5 8];

xtick_v = [x_limit(1),x_limit(1)+diff(x_limit)/4,...
           x_limit(1)+diff(x_limit)/2,...
           x_limit(1)+3*diff(x_limit)/4, x_limit(2)];
xticks(xtick_v);
xtick_v = xtick_v./1000;
xticklabels({num2str(xtick_v(1)),...
    num2str(xtick_v(2)),...
    num2str(xtick_v(3)),...
    num2str(xtick_v(4)),...
    num2str(xtick_v(5))});

ytick_v = [y_limit(1),y_limit(1)+diff(y_limit)/4,...
           y_limit(1)+diff(y_limit)/2,...
           y_limit(1)+3*diff(y_limit)/4,y_limit(2)];
yticks(ytick_v);
ytick_v = ytick_v./1000; %km
yticklabels({num2str(ytick_v(1)),...
    num2str(ytick_v(2)),...
    num2str(ytick_v(3)),...
    num2str(ytick_v(4)),...
    num2str(ytick_v(5))});

ztick_v = [z_limit(1),z_limit(1)+diff(z_limit)/3,...
    z_limit(1)+2*diff(z_limit)/3,z_limit(2)];
zticks(ztick_v);
ztick_v = ztick_v./1000; %km
zticklabels({num2str(ztick_v(1)),...
    sprintf('%0.1f',ztick_v(2)),...
    sprintf('%0.1f',ztick_v(3)),...
    num2str(ztick_v(4))});

xlabel('East (km)');
ylabel('North (km)');
zlabel('Altitude (km)');
title(sprintf('%s',datestr(obj.time)));

cellAll = [];
for plumeIdx=1:length(obj.plumes)
    plume = obj.plumes(plumeIdx);
    [vent_enu(1),vent_enu(2),vent_enu(3)] = ...
        plume.vent.enu(obj.areaLoc);
    plume_center = [vent_enu(1)-plume.res(1)/2,...
        vent_enu(2)-plume.res(2)/2,...
        vent_enu(3)];
    
    for rcvrIdx=1:length(obj.receivers)
        %Signal Paths
        for pathIdx=1:length(obj.signalPaths{rcvrIdx}(1,:))
            %Intersections        
            ray = obj.signalPaths{rcvrIdx}(pathIdx);              
            [isect_m, cell_m, line_v] = obj.signalIntersection(ray,plume);
            
%             if(size(isect_m)>0)    
% %                 plot3(isect_m(:,1),isect_m(:,2),isect_m(:,3),'ro',...
% %                     'MarkerSize',2,'markerfacecolor','r');
%                   %Plot intersect lines
%                 for ii=1:length(isect_m(:,1))-1
%                     line(isect_m(ii:ii+1,1),...
%                         isect_m(ii:ii+1,2),...
%                         isect_m(ii:ii+1,3),'Color','g','LineWidth',2);
%                 end
%             end
            cellAll = [cellAll;cell_m];
        end
    end
    
    plumePS = plume.getPixSize();
    x_v = linspace(-1*plume.size(1)/2+plume.res(1)/2,...
                   plume.size(1)/2-plume.res(1)/2,...
                   plumePS(1));
    y_v = linspace(-1*plume.size(2)/2+plume.res(2)/2,...
                   plume.size(2)/2-plume.res(2)/2,...
                   plumePS(2));
    z_v = linspace(0,plume.size(3)-plume.res(3),plumePS(3));
    
    
    for uu = 1:plumePS(3)
        for nn = 1:plumePS(2)
            cube_v = [];
            face_v = [];
            for ee = 1:plumePS(1)
                if(sum(plume.density(ee,nn,uu,:))>0)
                    cube = ...
                      [ver(:,1)*plume.res(1)+x_v(ee)+plume_center(1),...
                      ver(:,2)*plume.res(2)+y_v(nn)+plume_center(2),...
                      ver(:,3)*plume.res(3)+z_v(uu)+plume_center(3)];
                    

                    isInCell = sum(sum((cellAll(:,1)==ee & ...
                                        cellAll(:,2)==nn & ...
                                        cellAll(:,3)==uu)));

                    csize = size(cube_v);
                    if(isInCell>0)  
                        if(csize>0)
                        cMin = [min(cube_v(:,1)),...
                                min(cube_v(:,2)),...
                                min(cube_v(:,3))];
                            
                        cMax = [max(cube_v(:,1)),...
                                max(cube_v(:,2)),...
                                max(cube_v(:,3))];
                            
                        cubeComb = ...
                            [cMax(1),cMax(2),cMin(3);...
                             cMin(1),cMax(2),cMin(3);...
                             cMin(1),cMax(2),cMax(3);...
                             cMax(1),cMax(2),cMax(3);...
                             cMin(1),cMin(2),cMax(3);...
                             cMax(1),cMin(2),cMax(3);...
                             cMax(1),cMin(2),cMin(3);...
                             cMin(1),cMin(2),cMin(3)];
                         
%                         pt1 = patch('Faces',face_v,'Vertices',cube_v);
                        pt1 = patch('Faces',fac,'Vertices',cubeComb,'edgecolor','none');
%                         set(pt1,'edgecolor','none');
                        alpha(pt1,.1);
                        set(pt1,'FaceColor',[.4 .4 .4]);
                        end
                        pt = patch('Faces',fac,'Vertices',cube,'edgecolor','none');
%                         set(pt,'edgecolor','none');
%                         alpha(pt,.9);
                        set(pt,'FaceColor',[0 0 1]);
                    else    
                        face_v = [face_v;fac+csize(1)];
                        cube_v = [cube_v;cube];
%                         pt = patch('Faces',fac,'Vertices',cube);
%                         set(pt,'edgecolor','none');
%                         alpha(pt,.1);
%                         set(pt,'FaceColor',[.4 .4 .4]);
                    end
                end
                if(~isempty(cube_v))
                    cMin = [min(cube_v(:,1)),...
                        min(cube_v(:,2)),...
                        min(cube_v(:,3))];
                    
                    cMax = [max(cube_v(:,1)),...
                        max(cube_v(:,2)),...
                        max(cube_v(:,3))];
                    
                    cubeComb = ...
                        [cMax(1),cMax(2),cMin(3);...
                        cMin(1),cMax(2),cMin(3);...
                        cMin(1),cMax(2),cMax(3);...
                        cMax(1),cMax(2),cMax(3);...
                        cMin(1),cMin(2),cMax(3);...
                        cMax(1),cMin(2),cMax(3);...
                        cMax(1),cMin(2),cMin(3);...
                        cMin(1),cMin(2),cMin(3)];
                    
                    %                         pt1 = patch('Faces',face_v,'Vertices',cube_v);
                    pt1 = patch('Faces',fac,'Vertices',cubeComb);
                    set(pt1,'edgecolor','none');
                    alpha(pt1,.1);
                    set(pt1,'FaceColor',[.4 .4 .4]);
                end
            end
        end
        
    end
end


% view(45,30);
% xlim([-2000,4000]);
% ylim([-2000,4000]);
% view(0,90);

end

