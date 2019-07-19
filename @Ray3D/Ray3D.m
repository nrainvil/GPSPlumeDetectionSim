classdef Ray3D < handle
    %RAY3D Ray Object
    %   3D representation of a ray for ray tracing simulation
    
    properties
        origin = [0,0,0];   % (pixels) Origin of ray
        dir = [1,0,0];      % (unit vector) direction of ray
    end
    
    methods
        function obj=Ray3D(originIn,dirIn)
            % Ray3D Constructor for Ray3D object
            % Input:    origin (pixels) (optional)
            %           dir (unit vector) (optional)
            if(nargin>0)
                obj.origin = originIn;
            end
            if(nargin>1)
                obj.setDir(dirIn);
            end
        end
        
        function setDir(obj,dirIn)
            % setDir set the ray direction
            % Input:    dir (unit vector)
            obj.dir = dirIn./norm(dirIn);
        end
        
        function setAzEl(obj,azIn,elIn)
            % setAzEl set the ray direction by Azimuth and Elevation
            % Input:    azIn (deg) Azimuth, clockwise from North
            %           elIn (deg) Elevation, from local horizon
            
%             xR = sind(azIn).*cosd(elIn);
%             yR = cosd(azIn).*cosd(elIn);
%             zR = sind(elIn);
            
            azGeo = -1*azIn*pi()/180+pi()/2;
            [xR,yR,zR] = sph2cart(azGeo,elIn.*pi()/180,1);
            obj.setDir([xR,yR,zR]);
        end
        
        function [coords] = getXYZ(obj,r)
            % getXYZ find coordinates for ray path at a distance of r
            % Input:    r (pixels)
            % Output:   coords (pixels)
            coords = obj.origin + r.*obj.dir;
        end
        
        function [az,el] = getAzEl(obj)
            % getAzEl Return Az/El representation of ray direction
            % Output:   az (deg) Azimuth
            %           el (deg) Elevation
            [azR,elR,~] = cart2sph(obj.dir(1),obj.dir(2),obj.dir(3));
            el = elR.*180/pi();
            az = azR.*180/pi();
            if(az<0)
                az=az+360;
            end
        end
    end
    
end

