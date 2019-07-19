classdef Location < matlab.mixin.Copyable
    %LOCATION - ECEF Location
    
    properties
        ecef = [0 0 0];                      %ECEF X,Y,Z (m) 
        spheroid=wgs84Ellipsoid('meters');  %Default to wgs84
    end
    
    methods
        function obj=Location(lat,lon,alt)
            % Location  Constructor for Location Class
            % Input:    Lat  (Deg)
            %           Long (Deg)
            %           Alt  (m)
            if nargin >= 3
                obj.setLoc(lat,lon,alt);
            end
        end
        
        function setLoc(obj,lat,lon,alt)
            % setLoc  Set location
            % Input:    Lat  (Deg)
            %           Long (Deg)
            %           Alt  (m)
            [obj.ecef(1), obj.ecef(2), obj.ecef(3)] = ...
                geodetic2ecef(obj.spheroid,lat,lon,alt,'degrees');  
        end
        
        function setENU(obj,east,north,alt,centerLocation)
            % setLoc  Set location
            % Input:    east (m)
            %           north (m)
            %           up  (m)
            [clat,clon,calt] = centerLocation.lla();
            [obj.ecef(1), obj.ecef(2), obj.ecef(3)] = ...
                enu2ecef(east,north,alt,clat,clon,calt,obj.spheroid);
        end
        
        function [lat,lon,alt] = lla(obj)
            % lla Find Geodetic Coordinates from ECEF
            %   Output: Lat (Deg)
            %           Long (Deg)
            %           Alt (m)
            [lat,lon,alt] = ...
                ecef2geodetic(obj.spheroid,...
                                obj.ecef(1),obj.ecef(2),obj.ecef(3));
        end
        
        function [east,north,up] = enu(obj,centerLocation)
            % enu Find East,North,Up
            % Input:    xyz (m) - ECEF Center point
            % Output:   East (m)
            %           North (m)
            %           Up (m)
            [clat,clon,calt] = centerLocation.lla();
            [east,north,up] = ecef2enu(obj.ecef(1),obj.ecef(2),...
                                        obj.ecef(3),clat,clon,calt,...
                                        obj.spheroid,'degrees');
        end
        
        function shiftENU(obj,enuIn,centerLocation)
            % shiftENU Adjust location by east, north, up in m
            % Input:    enuIn (1x3 m)
            %           centerLocation (Location) - ECEF Center Point
            [lat,lon,alt] = centerLocation.lla();
            [enuOld(1),enuOld(2),enuOld(3)] = obj.enu(centerLocation);
            enuNew = enuIn + enuOld;
            [obj.ecef(1),obj.ecef(2),obj.ecef(3)] = ...
                            enu2ecef(enuNew(1),enuNew(2),enuNew(3),...
                                     lat,lon,alt,obj.spheroid,'degrees');
        end
    end
    
end

