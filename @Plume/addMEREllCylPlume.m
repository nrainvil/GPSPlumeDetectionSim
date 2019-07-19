function [refAlt] = addMEREllCylPlume(obj,saA,saB,dirA,topAlt,MER,Vvent,simLoc)
% addConicalPlume Create conical shaped plume
% Input:        saA (m) - semi-axis A
%               saB (m) - semi-axis B
%               dirA (2x1) - unit vector in direction of A axis
%               topAlt (m) - Altitude of top of plume
%               MER (kg/s) - Mass Eruption Rate
%               Vvent (m/s) - Velocity at vent
%               simLoc (Location) - Simulation Location

type = 3;
plumePS = obj.getPixSize();
g = 9.8; % m/s

tcnt = length(obj.diameter);
ps = obj.getPixSize();
obj.density = zeros(ps(1),ps(2),ps(3),tcnt);
[~,~,baseAlt] = obj.vent.enu(simLoc);
height = topAlt-baseAlt;
ellipseArea = pi()*saA*saB;
bVol = obj.res(1)*obj.res(2)*obj.res(3); % m^3


[~,~,vent_alt] = obj.vent.lla;
x_v = linspace(-1*obj.size(1)/2,obj.size(1)/2,plumePS(1)+1);
y_v = linspace(-1*obj.size(2)/2,obj.size(2)/2,plumePS(2)+1);

alt_v = linspace(vent_alt,topAlt,plumePS(3));
sliceHeight = alt_v(2)-alt_v(1);
blockCnt = 0;
massTot = 0;
for alt_idx=1:plumePS(3)
    alt = alt_v(alt_idx);
    refAlt = alt - vent_alt;
    %Find DRE volume based on MER
    if(alt_idx==1)
        Vstart = Vvent;
    end
    Vend = sqrt(2*(-g)*obj.res(3)+Vstart^2);
    %Check if Plume has reached max height
    if(~imag(Vend)>0)
        Vavg = (Vstart+Vend)/2; % m/s
        altDensity = (2*MER/Vavg)/ellipseArea; %kg/m^3
        DREMass = (MER/Vavg)*sliceHeight*1e3; %g
        DREMass = DREMass*2; %Rising and Falling ash
        %Weibull distribution, number of particles per size for
        %entire plume
        [sW, nmax] = psdSWMass(DREMass,obj.k,obj.uN,obj.diameter);
%         fprintf('Alt: %0.2f Density: %2.2e (kg/m^3) NMax: %2.2e (#) ',...
%             alt,altDensity,nmax);
        if(alt > baseAlt && refAlt < topAlt)
            for x_idx=1:length(x_v)
                for y_idx=1:length(y_v)
                    point = [x_v(x_idx),y_v(y_idx)];
                    if(inEllipse(point,saA,saB,dirA))
                        %Density g/m^3
                        blockCnt = blockCnt+1;
                        obj.density(x_idx,y_idx,alt_idx,:) = sW.*...
                            (4/3).*pi().*((obj.diameter/2).^3).*...
                            obj.particleDensity; % g/(Altitude Slice)
                    end
                end
            end
        end
    else
        fprintf('Max Plume Height: %0.2f ',refAlt);
        fprintf('Time airborne: %0.2f (s)\n',4*refAlt/(Vvent));
        break;
    end
    Vstart = Vend;
    obj.density(:,:,alt_idx,:) = ...
        obj.density(:,:,alt_idx,:)./blockCnt; % g/block
    massSlice = sum(obj.density(floor(plumePS(1)/2),...
                            floor(plumePS(2)/2),...
                            alt_idx,:))*blockCnt; % g
    massTot = massTot + massSlice;
%     fprintf('Calc Density: %2.2e (kg/m^3)',1e-3*massSlice/(blockCnt*bVol));
%     fprintf('Block Count: %d\n',blockCnt);
    blockCnt = 0;
end


obj.density = (obj.density)./bVol; % g/m^3
fprintf('Estimated Plume Mass: %0.2e (kg) ',MER*4*refAlt/(Vvent));
fprintf('Calculated Plume Mass: %0.2e (kg)\n',massTot*1e-3);

end

