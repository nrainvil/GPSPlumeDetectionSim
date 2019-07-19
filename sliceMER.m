classdef sliceMER < handle
    %SLICEMER Slice for Mass Eruption Rate sim
    
    properties
        diameter = 30;          % m
        height = 0;             % m
        enu = 3000;             % (1x3 m)
        mass = 100;             % g
        velocity = [0,0,10];    % m/s
    end
    
    methods
        function obj=sliceMER(diaIn,heightIn,enuIn,massIn,velIn)
           % Constructor for  sliceMER class
           % Input:     diaIn (m) diameter of slice
           %            heightIn (m) height of slice
           %            enuIn (m) position of slice
           %            massIn (g) mass of slice
           %            velocity (1x3 m/s) velocity of slice in ENU
           obj.diameter = diaIn;
           obj.height = heightIn;
           obj.enu = enuIn;
           obj.mass = massIn;
           obj.velocity = velIn;
        end
    end
    
end

