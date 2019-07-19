classdef VaprSim < handle
    %VAPRSIM Simulation Environment
    %   Top level Ash Plume Attenuation Simulation Environment
    
    properties
        areaLoc = Location(37.734,15.004,0);  % Center Location
        areaSize = [8000,8000,8000];     % Size (m)
        areaRes = [100,100,100];            % Resolution (m)
        time = datetime(2018,1,1,1,0,0);    % (datetime) Current time 
        receivers = [];                     % (VaprRcvr Cell Array)
        plumes = [];                        % (Plume Cell Array)
        vents = [];                         % (Location Cell Array)
        gps = GPSCon();                     % (GPSCon) GPS Constellation
        signalPaths = {};                   % (ray3D Cell Array) 
        signalPRNs = {};                    % (# Cell Array) 
        signalAtten = {};
        signalLen = {};
        demX = [];
        demY = [];
        demZ = [];
        freq = 1.57542e9;                  %Hz L1
%         freq = 1.2276e9;                  %Hz L2
        minEle = 5;                         % (deg)
        plotLim = [-2000,4000];
        plotView = [45,30];
        mie = 1;                            % Set to 1 to enable mie x-s 
                                            % 2 for charged, 3 for csv
        eta = 0;                       % Surface charge
    end
    
    properties (Access = private)   
    end
    
    properties (Constant)
        u0 = 4e-7*pi();                     %Vacuum permeability
        e0 = 8.854187817e-12;               %Vacuum permittivity
        Np2dB = 10*log10(exp(1));       %Np/m to dB/m conversion, 4.15
    end
    
    methods
        [fig] = plotSim(obj,hideFig);
        [xB,yB,zB] = signalBoundary(obj, signalIn);
        [isect_m, cell_m, length_v] = ...
                                    signalIntersection(obj,rayIn, plumeIn);
        voldoradG3(obj);
    end
    
    methods
        function obj=VaprSim(timeIn,areaSizeIn,areaResIn)
            % VaprSim Constructor for VaprSim Object
            % Input:    areaSizeIn [x, y, z] (m) (optional)
            %           areaResIn [x, y, z] (m)  (optional)
            if nargin >=1
                obj.time = timeIn;
            end
            if nargin >=2
                obj.areaSize = areaSizeIn;
            end
            if nargin >=3
                obj.areaRes = areaResIn;
            end
            obj.setTime(obj.time);
        end
        
        function addRec(obj,recObj)
            % addRec  Add Receiver to Simulation
            % Input: recObj (VaprRcvr)
            obj.receivers = [obj.receivers, recObj];
        end
        
        function addVent(obj,locObj)
            % addVent  Add Vent to Simulation
            % Input: locObj (Location)
            obj.vents = [obj.vents, locObj];
        end
        
        function addPlume(obj,plumeObj)
            % addPlume  Add Plume to Simulation
            % Input: plumeObj (Plume)
            obj.plumes = [obj.plumes, plumeObj];
        end
        
        function [size_v] = getPixSize(obj)
            % getPixScale Return the scaling factor for a pixel based
            % ray tracing simulation
            % Output:   sX (m/pixel)
            %           sY (m/pixel)
            %           sZ (m/pixel)
            size_v = obj.areaSize./obj.areaRes; 
        end
        
        function setTime(obj,timeIn)
           % setTime  Change time of simulation 
           % Input: timeIn (datetime)
           obj.time=timeIn;
           obj.gps.findOrbits(timeIn);
           obj.findSignalPaths(obj.minEle);
           %obj.findSignalAtten();
        end
        
        function findSignalPaths(obj,minEle)
           % findSignalPaths  Find the receiver to satellite paths 
           % Input: minEle (Deg) Minimum Elevation Angle
           obj.signalPaths = [];
           obj.signalPRNs = [];
           for rcvrIdx=1:length(obj.receivers)
               rcvr = obj.receivers(rcvrIdx);
               [origin(1),origin(2),origin(3)] = rcvr.enu(obj.areaLoc);
               [az_v,el_v,sat_v] = obj.gps.azelGTEle(rcvr,minEle);
               
               for ii=1:length(sat_v)
                   %Create ray with origin at ENU location
                    ray = Ray3D(origin./obj.areaRes,[1,0,0]);
                    %Set ray direction based on azimuth degrees from North
                    %and elevation, degrees from local horizon.
                    ray.setAzEl(az_v(ii),el_v(ii));
                    obj.signalPaths{rcvrIdx}(ii)=ray;
                    obj.signalPRNs{rcvrIdx}(ii) = sat_v(ii);
               end
           end
           
        end
        
        function [pathMass,pathLen,atn] = findSignalAtten(obj,rayIn,debug)
           % findSignalAtten  Find the attenuation for each signal
           % Output:        pathMass (Mx1 kg/m^2 ?)
           %                atn (dB)
           
           if(obj.mie>2)
            fid = fopen('vsr_sim.csv','a');
           end
           
           if (nargin<3)
               debug=0;
           end
           
           for plumeIdx=1:length(obj.plumes)
               %Find Plume plane intersections along ray path
               [~, cell_m, length_v] = ...
                        obj.signalIntersection(rayIn,obj.plumes(plumeIdx));
               dens_m = obj.plumes(plumeIdx).density; % g/m^3
%                ndens_m = obj.plumes(plumeIdx).nDensity; % #/m^3
               inP_m = dens_m;
               inP_m(inP_m>0) = 1;
               
               volTot = 0;
               massTot = 0;
               atn = 0;
               atn_a = 0;
               pathMass = 0;
               pathLen = 0;
               %Iterate over ash types
               for typeIdx=1:length(obj.plumes(plumeIdx).diameter)
                   densIdx_v = sub2ind(size(dens_m),...
                       cell_m(:,1),cell_m(:,2),cell_m(:,3),...
                       typeIdx*ones(size(cell_m(:,1))));
                   mass_v = dens_m(densIdx_v).*length_v; % g/m^2
                   pathMass = pathMass + sum(mass_v);
%                    pathLen = ...
%                             max([pathLen,sum(inP_m(densIdx_v).*length_v)]);
                   inPlumeMask = inP_m(densIdx_v)>0;
%                    if(sum(inPlumeMask)>0&&typeIdx==1)
%                         for kk=1:length(cell_m(:,1))
%                             if(inPlumeMask(kk)>0)
%                             fprintf('%0.2f %0.2f %0.2f|%0.2f\n',...
%                                         cell_m(kk,1),cell_m(kk,2),...
%                                         cell_m(kk,3),length_v(kk));
%                             end
%                         end
%                    end
                   pathLen = max([pathLen,sum(length_v(inPlumeMask))]); % m
                   
                   %Calculate Extinction Coeff
                   rho = obj.plumes(plumeIdx).particleDensity(typeIdx); % g/m^3
                   dia = obj.plumes(plumeIdx).diameter(typeIdx); % m
                   eSct = obj.plumes(plumeIdx).permittivity(typeIdx);
                   
                   if(obj.mie==0)
                    c_a = cross_abs(dia/2, obj.freq, eSct);
                    c_s = cross_sct(dia/2, obj.freq, eSct);
                   elseif(obj.mie==1)
                    c_a = mie_cross_abs(dia/2, obj.freq, eSct);
                    c_s = mie_cross_sct(dia/2, obj.freq, eSct);
                   elseif(obj.mie==2)  
                    c_a = chrg_cross_abs(dia/2, obj.freq, eSct, obj.eta);
                    c_s = chrg_cross_sct(dia/2, obj.freq, eSct, obj.eta);
                   else
                    fprintf(fid,'%0.8f,%0.8f,%0.8f,%0.8f\n',dia,obj.freq,real(eSct),imag(eSct));
                    c_a = cross_abs(dia/2, obj.freq, eSct);
                    c_s = cross_sct(dia/2, obj.freq, eSct);
                   end
                   pvol = (4/3)*pi()*(dia/2).^3; % m^3
                   pmass = rho.*pvol; % g
                   ndens_v = dens_m(densIdx_v)./pmass; % #/m^3
%                    ndens_v2 = ndens_m(densIdx_v); % #/m^3
                   if(debug>=1)
                   volTot = volTot + max(ndens_v);
                   massTot = massTot + max(dens_m(densIdx_v));
                   end
                   if(debug>=2)
                       if(max(ndens_v)>0)
                   fprintf('%0.4e(#/m^3) %0.2e (g/m^3) ',...
                            max(ndens_v),max(dens_m(densIdx_v)));
                   fprintf('d: %0.2f (mm) c_a: %0.4e (m) c_s: %0.4e (m) l: %0.4e (m)\n',...
                            dia*1e3,c_a,c_s,sum(length_v));
                       end
                   end
                   
                   %Calculate Extinction Coefficient
                   c_e = c_a + c_s;
                   ke_v = ndens_v.*(c_e); %Np/m 
                   ke_v = obj.Np2dB.*ke_v;  %dB/m
                   ke_v(isnan(ke_v)) = 0;
                   %Find attenuation
                   atn = atn + sum(ke_v.*length_v);
                   
                   ke_a_v = obj.Np2dB.*ndens_v.*(c_a);
                   ke_a_v(isnan(ke_a_v)) = 0;
                   atn_a = atn_a + sum(ke_a_v.*length_v);
               end
%                atn = (1/20)*log(10)*atn; %N to dB
%                atn_a = (1/20)*log(10)*atn_a; %N to dB
               if(debug>=1 && atn>0)
                   fprintf(' Absorption: %0.4f (dB) ',atn_a);
                   fprintf('Scattering: %0.4f (dB) ',atn-atn_a);
                   fprintf('Extinction: %0.4f (dB)\n',atn);
                   fprintf(' Total Volume Density: %0.2f (#/m^3) ',volTot);
                   fprintf('Mass Density: %0.2f (g/m^3)\n',massTot);
               end
           end
           if(obj.mie>2)
               fclose(fid);
           end
        end
                
        function [isInSim]=inSim(obj,coords_v)
            % inSim Find if point is in simulation area
            % Input:    coords_v (1x3 m) ENU
            % Output:   isInSim (bool) 1 if point is in sim
            isInSim = ~((abs(coords_v(1)) > obj.areaSize(1)/2) || ...
                     (abs(coords_v(2)) > obj.areaSize(2)/2) || ...
                     (coords_v(3) > obj.areaSize(3)) || (coords_v(3) < 0));
            if(isnan(coords_v))
                isInSim = 0;
            end
        end
        
        function loadDEM(obj,filePath)
            % loadDEM  Load DEM file
            % Input: filePath (string) 
            fid = fopen(filePath,'r');
            csv = textscan(fid,'%f%f%f','Delimiter',',','Headerlines',1);
            dem_lla=[csv{1},csv{2},csv{3}];
            [area_lat,area_lon,area_alt] = obj.areaLoc.lla;
            [demFull(:,1),demFull(:,2),demFull(:,3)] = ...
                geodetic2enu(dem_lla(:,1),dem_lla(:,2),dem_lla(:,3),...
                             area_lat,area_lon,area_alt,...
                             obj.areaLoc.spheroid,'degrees');
            
            demSize(1) = max(demFull(:,1)) - min(demFull(:,1));            
            demSize(2) = max(demFull(:,2)) - min(demFull(:,2));
            
            %Select DEM inside area size             
                         
            %Resample to resolution             
            [obj.demX,obj.demY] = ...
                 meshgrid(linspace(min(demFull(:,1)),...
                                   max(demFull(:,1)),...
                                   demSize(1)./obj.areaRes(1)),...
                          linspace(min(demFull(:,2)),...
                                   max(demFull(:,2)),...
                                   demSize(2)./obj.areaRes(2)));
            obj.demZ = griddata(demFull(:,1),demFull(:,2),demFull(:,3),...
                                obj.demX,obj.demY);                   
            fclose(fid);
        end
      
        function altOut = findAlt(obj,eIn,nIn)
        % findAlt Find Altitude from DEM
        % Input: eIn (m) East
        %        nIn (m) North
        % Output: altOut (m) Altitude
            [~,minELoc] = min(abs(obj.demX(1,:)-eIn));
            [~,minNLoc] = min(abs(obj.demY(:,1)-nIn));
            altOut = obj.demZ(minNLoc,minELoc);
        end
    end
    
end

