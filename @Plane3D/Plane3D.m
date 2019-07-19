classdef Plane3D < handle
    %PLANE3D Plane object
    %   3D representation of a plane for ray tracing simulation
    
    properties
        point = [0,0,0];    % (pixels) point on the plane
        normal = [0,0,1];   % (unit vector) vector normal to plane
    end
    
    methods
        function obj=Plane3D(pointIn,normalIn)
        % Plane3D Contructor for Plane3D object
        % Input:    pointIn (pixels) (optional)
        %           normalIn (unit vector) (optional)
            if(nargin>0)
                obj.point = pointIn;
            end
            if(nargin>1)
                obj.setNormal(normalIn);
            end
        end
        
        function setNormal(obj,normalIn)
        % setNormal set the normal vector of the plane
        % Input:    normalIn (unit vector)
            obj.normal = normalIn./norm(normalIn);
        end
        
        function [dist]=rayIntersect(obj,rayIn)
        % rayIntersect Find intersection of ray with plane
        % Input:    rayIn (Ray3D)
        % Output:   dist (pixels) distance on ray to intersection
        
        rayDiv = dot(rayIn.dir,obj.normal);
        dist = dot(obj.point-rayIn.origin,obj.normal)./rayDiv;
        
        if((dist < 0) || (abs(rayDiv) < 1e-6))
            dist = NaN;
        end
        end
    end
    
end

