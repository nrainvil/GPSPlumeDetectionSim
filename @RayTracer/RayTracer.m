classdef RayTracer < handle
    %RAYTRACER Ray Tracer
    %   Ray Tracing for emag simulations
    
    properties
        simSize = [120,120,120];     % Pixels
    end
    
    methods
        function obj=RayTracer()
        % RayTracer Constructor for Ray Tracer Object
        end
    end
    
    methods(Static)
        function [r] = intersectPlane(rayIn,planeIn)
        % intersectPlane Calculate the distance from the ray origin to the
        % intersection with the plane
        % Input:    rayIn (Ray3D object)
        %           planeIn (Plane3D object)
        % Output:   r (pixels) Distance to the plane
            projRP = dot(rayIn.dir,planeIn.normal);
            projRP(abs(projRP) < 1e-6) = NaN;
            r = dot((planeIn.point-rayIn.origin),planeIn.normal)./projRP;
            r(r<0) = NaN;
        end
        
    end
    
end

