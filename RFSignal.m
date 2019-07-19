classdef RFSignal < Ray3D
    %RFSignal RF Signal Path
    %   3D representation of an RF Signal
    
    properties
        scale = 1;     % (m/pixel)
    end
    
    methods
        function obj=RFSignal(scaleIn, originENU, azIn, elIn)
            %RFSignal Constructor for RF Signal
            % Input:    scaleIn (m/pixel) scaling factor for simulation
            %           originENU ([m,m,m]) XYZ in enu coordinates
            %           azIn (deg) azimuth of signal path
            %           elIn (deg) elevation of signal Path
            
            obj@Ray3D();
            
            if(nargin>0)
                obj.scale=scaleIn;
            end
            
            if(nargin>2)
               zR =  originENU(3);
               r_v = zR./tand(elIn);
               xR = r_v.*sind(azIn)+originENU(1);
               yR = r_v.*cosd(azIn)+originENU(2);
               
               obj.origin = originENU./obj.scale;
               obj.setDir(xR,yR,zR);
            end
        end
        
        function [coords] = getXYZ(obj,r)
            % getXYZ find coordinates for signal path at a distance of r
            % Input:    r (m)
            % Output:   coords (m) ENU
            coords = obj.scale.*(getXYZ@Ray3D(obj,r./obj.scale));
        end
    end
end

