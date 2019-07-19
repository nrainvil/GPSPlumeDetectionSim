classdef Plume < handle
    %PLUME Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vent = Location();          % Center location
        size = [6000,6000,6000];    % (m)
        res = [100,100,100];        % (m)
        density                     % (MxNxVxTypes) density (g/m^3)
%         nDensity                    % (MxNxVxTypes) density (#/m^3)
        particleDensity             % (Typesx1) density (g/m^3)
        diameter                    % (Typesx1) (mm)
        permittivity                % (Typesx1) permitttivity 
        k = 2.3;                      % PSD Shape Factor
        uN = .0164;                 % (m) PSD Particle Size Mode
        slices_v = {};              % (Timex1 sliceMER) For MER sim
        relMass = [];               % Relative Mass of each particle
    end
        
    methods
        [ext_coeff] = getExtCoeff(obj,freqIn,cellX,cellY,cellZ);
        addConicalPlume(obj,baseDia,topDia,topAlt,DREVolume,simLoc);
        addEllipticalCylPlume(obj,saA,saB,dirA,topAlt,DREVolume,simLoc);
        [refAlt] = addMEREllCylPlume(obj,saA,saB,dirA,topAlt,MER,Vvent,simLoc);
    end
    
    methods
        function obj=Plume(ventIn,sizeIn,resIn)
            % Plume  Constructor for Plume Class
            % Input:    vent (Location) (optional)
            %           sizeIn [x,y,z] (m) (optional)
            %           res [x,y,z] (m) (optional)
            if nargin >= 1
                obj.vent = ventIn;
            end
            if nargin >=2
                obj.size = sizeIn;
            end
            if nargin >=3
                obj.res = resIn;
            end
            
            obj.size = obj.res.*obj.getPixSize();
        end
        
        function initParticleTypes(obj,typeDia_v,typeE_v,typeRho_v)
            % initParticleTypes Initialize Particle Types
            % Input:        typeDia_v = particle diamter (m)
            %               typeE_v = complex particle permittivity
            %               typeRho_v = particle density (kg/m^3)
            ps = obj.getPixSize();
            tcnt = length(typeDia_v);
            
            obj.particleDensity = typeRho_v;
            obj.permittivity = typeE_v;
            obj.diameter = typeDia_v;
            obj.density = zeros(ps(1),ps(2),ps(3),tcnt);
%             obj.nDensity = zeros(ps(1),ps(2),ps(3),tcnt);
        end
        
        function [pixSize]=getPixSize(obj)
           % getPixSize Return the size of the plume in pixels (m/m)
           % Output:    pixSize (1x3) enu
           pixSize = ceil(obj.size./obj.res);
        end
        
        function [isInPlume]=inPlume(obj,coords_v,sim)
            % inPlume Find if point is in Plume area
            % Input:    coords_v (1x3 m) ENU
            % Output:   isInPlume (bool) 1 if point is in plume
            [center_v(1),center_v(2),center_v(3)] = ...
                                    obj.vent.enu(sim.areaLoc);
            coordsC_v = coords_v - center_v;
            isInPlume = ~((abs(coordsC_v(1)) > obj.size(1)/2) || ...
                     (abs(coordsC_v(2)) > obj.size(2)/2) || ...
                     (coordsC_v(3) > obj.size(3)) || (coordsC_v(3) < 0));
            if(isnan(coordsC_v))
                isInPlume = 0;
            end
%             if(isInPlume)
%                fprintf('IN: %0.2f,%0.2f,%0.2f\n',...
%                    coords_v(1),coords_v(2),coords_v(3));
%             else
%                 fprintf('OUT:  %0.2f,%0.2f,%0.2f -- %0.2f,%0.2f,%0.2f\n',...
%                    coords_v(1),coords_v(2),coords_v(3),...
%                    center_v(1),center_v(2),center_v(3));
%             end
        end
        
        function pslc=plot(obj,vslice)
%            f = figure('visible','off');
           vslice_v = linspace(1/(vslice+1),vslice/(vslice+1),vslice);
           pslc = slice(obj.density,obj.res(1)/2,obj.res(2)/2,...
                    obj.res(3)*vslice_v); 
                
           set(pslc(1),'LineWidth',.01);       
           set(pslc(2),'LineWidth',.01);
           for ii=3:length(pslc)
                set(pslc(ii),'edgecolor','none');
           end

           
        end
    end
    
end

