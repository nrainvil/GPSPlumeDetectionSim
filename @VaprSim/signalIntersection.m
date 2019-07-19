function [isect_m, cell_m, length_v] = ...
                                    signalIntersection(obj,rayIn, plumeIn)
% rayIntersection Calculate intersection of ray with Plume
% Input:    rayIn (ray3D)
%           plumeIn (Plume)
% Output:   isect_v (M+1x3) Intersection points
%           cell_m  (Mx3) Location of intersected cells
%           length_v (Mx1) Length of paths through cells
cellCnt = plumeIn.getPixSize();
[plumeOrigin(1), plumeOrigin(2), plumeOrigin(3),] = ...
    plumeIn.vent.enu(obj.areaLoc);
isect_m = [];
cell_m = [];
length_v = [];
pPlanes_v = [];

%Create Plume Planes
for ee = 1:(cellCnt(1)+1)
    eastLoc = (ee-((cellCnt(1)+1)/2)-.5).*plumeIn.res(1);
    loc_v = [eastLoc,0,0] + plumeOrigin;
    tPlane = Plane3D(loc_v./obj.areaRes,[1,0,0]); % east
    pPlanes_v = [pPlanes_v;tPlane];
end

for nn = 1:(cellCnt(2)+1)
    northLoc = (nn-((cellCnt(2)+1)/2)-.5).*plumeIn.res(2);
    loc_v = [0,northLoc,0] + plumeOrigin;
    tPlane = Plane3D(loc_v./obj.areaRes,[0,1,0]); % north
    pPlanes_v = [pPlanes_v;tPlane];
end

for uu =1:(cellCnt(3)+1)
    upLoc = (uu-1).*plumeIn.res(3);
    loc_v = [0,0,upLoc] + plumeOrigin;
    tPlane = Plane3D(loc_v./obj.areaRes,[0,0,1]); % up
    pPlanes_v = [pPlanes_v;tPlane];
end

%Check if intersection is in Plume area
for ii=1:length(pPlanes_v)
    tPlane = pPlanes_v(ii);
    distPix = tPlane.rayIntersect(rayIn);
    coords = rayIn.getXYZ(distPix).*obj.areaRes; %enu
    if(plumeIn.inPlume(coords,obj))
        isect_m = [isect_m; coords];
%         fprintf('%d:%0.2f,%0.2f,%0.2f\n',ii,coords(1),coords(2),coords(3));
    end
end

%Find intersected Cells
for ee = 1:cellCnt(1)
    for nn = 1:cellCnt(2)
        for uu = 1:cellCnt(3)
            %Create Faces
            cCenter =[(ee-((cellCnt(1)+1)/2)).*plumeIn.res(1),...
                         (nn-((cellCnt(2)+1)/2)).*plumeIn.res(2),...
                         (uu-.5).*plumeIn.res(3)]...
                         + plumeOrigin;
             face_c{1}(ee,nn,uu) = cCenter(1)+plumeIn.res(1)/2; % east
             face_c{2}(ee,nn,uu) = cCenter(1)-plumeIn.res(1)/2; % west
             face_c{3}(ee,nn,uu) = cCenter(2)+plumeIn.res(2)/2; % north
             face_c{4}(ee,nn,uu) = cCenter(2)-plumeIn.res(2)/2; % south
             face_c{5}(ee,nn,uu) = cCenter(3)+plumeIn.res(3)/2; % up
             face_c{6}(ee,nn,uu) = cCenter(3)-plumeIn.res(3)/2; % down
             %Check of intersection is in cell
        end
    end
end

round12 = @(x) round(x,10); %Round to 12 digits of precision
face_c = cellfun(round12,face_c,'UniformOutput',false);

plumeMask = zeros(cellCnt);
limTest = [2,3;
           3,1;
           1,2];
res_m = [];
if(~isempty(isect_m))
for ii=1:length(isect_m(:,1))
    isect = isect_m(ii,:);
    %Test Each Face
    for jj=1:6
        dirInd = ceil(jj/2);
        ind = find(...
            face_c{2*limTest(dirInd,1)-1}>isect(limTest(dirInd,1)) &...
            face_c{2*limTest(dirInd,1)}<isect(limTest(dirInd,1)) &...
            face_c{2*limTest(dirInd,2)-1}>isect(limTest(dirInd,2)) &...
            face_c{2*limTest(dirInd,2)}<isect(limTest(dirInd,2)) &...
            face_c{jj}==round12(isect(dirInd)));
        if(~isempty(ind))
            [e,n,u] = ind2sub(size(face_c{1}),ind);
%             fprintf('%d|%d: %d,%d,%d\n',ii,jj,e,n,u);
            res_m= [res_m;[ii,e,n,u]]; 
            plumeMask(e,n,u) = plumeMask(e,n,u) + 1;
        end
    end
end
end

%Find Length
[e_v,n_v,u_v] = ind2sub(size(plumeMask),find(plumeMask>1));
cell_m = [e_v,n_v,u_v];
if(~isempty(cell_m))
    length_v = zeros(length(cell_m(:,1)),1);
    for ii=1:length(cell_m(:,1))
        resCell = res_m(res_m(:,2)==cell_m(ii,1) &...
            res_m(:,3)==cell_m(ii,2) &...
            res_m(:,4)==cell_m(ii,3));
        length_v(ii,1) = norm(isect_m(resCell(1),:)-isect_m(resCell(2),:));
%         fprintf('Isect1: %0.2f %0.2f %0.2f',isect_m(resCell(1),1),...
%                     isect_m(resCell(1),2),isect_m(resCell(1),3));
%         fprintf(' Isect2: %0.2f %0.2f %0.2f',isect_m(resCell(2),1),...
%                     isect_m(resCell(2),2),isect_m(resCell(2),3));
%         fprintf(' Length: %0.2f\n',length_v(ii,1));
    end
end
end

