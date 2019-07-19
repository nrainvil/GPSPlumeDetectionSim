classdef GPSCon < handle
    %GPSCON GPS Constellation Model
    
    properties
        time = datetime(2018,1,1);  %Current Time 
        sat_v = [];                 %Satellite Number
        ecef_v = [];                %ECEF Position of satellite
        sp3 = [];
    end
    
    methods
        function obj=GPSCon(timeIn)
            % GPSCon  Constructor for GPSCon Class
            % Input:    timeIn (datetime) (optional)
            if nargin >=1
                obj.findOrbits(timeIn);
            end
        end
        
        function findOrbits(obj,time)
            % findOrbits  Set the location of all satellites to the current
            %             time.
            % Input:    time (datetime)
            
            maxConst = 1; % 1 for just GPS, 2 for GPS + GLONASS, 
            
            if(isempty(obj.sp3) || ...
                    (obj.time.Day ~= time.Day) || ...
                    (obj.time.Month ~=obj.time.Month) || ...
                    (obj.time.Year ~= obj.time.Year))
                
                obj.time = time;
                [obj.sp3, found_orbits] = ...
                    load_sp3_information(obj.time.Year,obj.time.Month,...
                    obj.time.Day);
%                 fprintf('Loading SP3 for %d/%d/%d\n',...
%                     obj.time.Month,obj.time.Day,obj.time.Year);
                if(~found_orbits)
                    fprintf('Orbits not found\n');
                    return;
                end
            else
                obj.time = time;
            end
            
            orbitHour = obj.time.Hour;
            [obj.sat_v, obj.ecef_v, ~] = ...
                do_orbits(obj.sp3,obj.time.Year,obj.time.Month,obj.time.Day,...
                            orbitHour,obj.time.Minute,maxConst);
            
        end
        
        function [az,el] = azelSat(obj,rcvr,sat_num)
            % azelSat  Return the azimuth and elevation for a satellite
            % Input:    rcvr (VaprRcvr) - Receiver
            %           sat_num - PRN Number
            % Output:   Az (Deg)
            %           El (Deg)
            
            sat_ecef = obj.ecef_v(obj.sat_v==sat_num,:);
            [lat,lon,~] = rcvr.lla;
            [az,el] = ecef2azelrange(sat_ecef',rcvr.ecef',lat,lon); 
        end
        
        function [az_v,el_v,sat_v] = azelGTEle(obj,rcvr,ele)
            % azelGTEle  Return the azimuth and elevation for all
            %            satellites above a reference elevation
            % Input:    rcvr (VaprRcvr) - Receiver
            %           ele (Deg) - Reference elevation
            % Output:   Az (Deg)
            %           El (Deg)
            %           sat (#) - PRN
            
            az_v = zeros(1,length(obj.sat_v));
            el_v = zeros(1,length(obj.sat_v));
            for prn_idx=1:length(obj.sat_v)
                prn = obj.sat_v(prn_idx);
                [az_v(prn_idx), el_v(prn_idx)] = obj.azelSat(rcvr,prn);
            end
            sat_v = obj.sat_v(el_v>=ele);
            az_v = az_v(el_v>=ele);
            el_v = el_v(el_v>=ele);
        end
    end
    
end