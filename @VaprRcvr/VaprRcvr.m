classdef VaprRcvr < Location
    %VAPRRCVR - VAPR Receiver Model
    
    properties
        antHeight = 1;  % (m) Antenna height
        name = 'NONE';  % Receiver name
    end
    
    methods
        function obj = VaprRcvr(lat,lon,alt,nameIn)
            % VaprRcvr  Constructor for VaprRcvr Class
            % Input:    Lat  (Deg)
            %           Long (Deg)
            %           Alt  (m)
            obj@Location();
            obj.setLoc(lat,lon,alt+obj.antHeight);
            if(nargin>3)
                obj.setName(nameIn);
            end
        end
        
        function setName(obj,nameIn)
            obj.name=nameIn;
        end
    end
    
end

