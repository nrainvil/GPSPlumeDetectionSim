function addEllipticalCylPlume(obj,saA,saB,dirA,topAlt,DREVolume,simLoc)
% addConicalPlume Create conical shaped plume
% Input:        saA (m) - semi-axis A
%               saB (m) - semi-axis B
%               dirA (2x1) - unit vector in direction of A axis
%               topAlt (m) - Altitude of top of plume
%               DREVolume (m^3) - Dense Rock Equivalent Volume
%               simLoc (Location) - Simulation Location

type = 3;
plumePS = obj.getPixSize();

[~,~,baseAlt] = obj.vent.enu(simLoc);
height = topAlt-baseAlt;
ellipseArea = pi()*saA*saB;
plumeVolume = ellipseArea*height; % m^3

meanDensity = (DREVolume*obj.particleDensity(type))/plumeVolume; % g/m^3
fprintf('Plume Volume: %0.2e m^3 ',plumeVolume);
fprintf('Mean Density: %0.2f (g/m^3) ',meanDensity);
fprintf('Plume Mass: %0.2e (kg)\n',plumeVolume*meanDensity*1e-3);

%Weibull distribution, number of particles per size for
%entire plume
[sW, nmax] = psdSWVol(DREVolume,obj.k,obj.uN,obj.diameter);
fprintf('NMax: %4.4e ',nmax);
[~,~,vent_alt] = obj.vent.lla;
x_v = linspace(-1*obj.size(1)/2,obj.size(1)/2,plumePS(1)+1);
y_v = linspace(-1*obj.size(2)/2,obj.size(2)/2,plumePS(2)+1);

alt_v = linspace(vent_alt,topAlt,plumePS(3));
blockCnt = 0;
for alt_idx=1:plumePS(3)
    alt = alt_v(alt_idx);
    refAlt = alt - vent_alt;
    if(alt > baseAlt && refAlt < topAlt)
        for x_idx=1:length(x_v)
            for y_idx=1:length(y_v)
                point = [x_v(x_idx),y_v(y_idx)];
                if(inEllipse(point,saA,saB,dirA))
                    %Density g/m^3
                    blockCnt = blockCnt+1;
                    obj.density(x_idx,y_idx,alt_idx,:) = sW.*...
                        (4/3).*pi().*((obj.diameter/2).^3).*...
                        obj.particleDensity; % kg/Plume
                end
            end
        end
    end
end

bVol = obj.res(1)*obj.res(2)*obj.res(3); % m^3
fprintf('Block Count: %d Plume Volume: %0.2e',blockCnt,blockCnt*bVol);
obj.density = obj.density./blockCnt; % g/block
obj.density = (obj.density)./bVol; % g/m^3

bDens = sum(obj.density(floor(plumePS(1)/2),...
                        floor(plumePS(2)/2),...
                        floor(plumePS(3)/2),:))*1e-3; %kg/m^3
                    
fprintf(' Plume Mass: %0.2e (kg)\n',blockCnt*bVol*bDens);

end

