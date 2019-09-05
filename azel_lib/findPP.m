function [edgeDist, deltaHeight, ppENU ] = ...
    findPP( vent, rcvr, azel, plumeRad, plumeHeight, center)
%FINDPP Find Pierce Point
%Input:     vent (Location)
%           rcvr (VaprRcvr)
%           azel ([deg deg]) GPS Azimuth
%           plume_rad (m) Radius of plume
%           center (Location)
%Output:    edgeDist (m) in East/North plane
%           deltaHeight (m) relative to vent height
%           ppENU (East,North,Up m) 

    azIn = mod(360+azel(:,1),360);
    elIn = mod(360+azel(:,2),360);
    
    ppENU = [0,0,0];
    min_h = 500; %m above vent
    elVent = [];
    deltaHeight = 0;
    edgeDist = 0;
    maxEle = 0;
   
    %Recevier Location
    [rcvrENU(1), rcvrENU(2), rcvrENU(3)] = rcvr.enu(center);
    %Vent Location
    [ventENU(1), ventENU(2), ventENU(3)] = vent.enu(center);
    [azLim,~] = findAzLim(ventENU, rcvrENU, plumeRad);
    azLim = gpsAz(azLim); %Azimuth Limits
    
    %Find Center
    if((azLim(2)-azLim(1))>180)
        azLim(2) = azLim(2)+360;
    end
    azCenter = azLim(1)+(azLim(2)-azLim(1))/2;
    
%     fprintf(' %s Az In: %0.2f Az Lim %0.2f|%0.2f|%0.2f ',...
%         rcvr.name,azIn,azLim(1),azCenter,azLim(2));
    
    %Elevation
    minVENU = ventENU(:) + [0;0;min_h];
    elVent = asind((minVENU(3)-rcvrENU(3))/norm(minVENU-rcvrENU));
    
    if(azIn > azLim(2) || azIn < azLim(1))
%         fprintf('Outside limits\n');
        return;
    end
    
    centerDist = abs(norm(ventENU(1:2)-rcvrENU(1:2)));
    
    %Edge of plume
    azHalfWidth = azCenter-azLim(1); %Relative center az
    azFlag = 0;
    az_ref = azHalfWidth-(azIn-azLim(1));
    if(az_ref<0)
        az_ref = azHalfWidth-(azLim(2)-azIn);
        azFlag = 1;
    end

    if(az_ref==0)
        edgeDist=centerDist-plumeRad;
    else
        th3 = 180-az_ref-asind(centerDist.*((plumeRad./sind(az_ref))^(-1)));
        th3(th3>90) = 180-2*az_ref-th3;
        edgeDist = (sind(th3)./sind(az_ref))*plumeRad;
    end
    
    plumeRelHeight = plumeHeight + (ventENU(3)-rcvrENU(3));
    maxEle = atand(plumeRelHeight/edgeDist);
%     fprintf(' Az Ref: %0.2f Max El:%0.2f \n',az_ref,maxEle);
    
    minEle = atand((ventENU(3)-rcvrENU(3))/centerDist);
    deltaHeight = edgeDist*tand(elIn)-ventENU(3)+rcvrENU(3);  
    if(elIn>maxEle | elIn < minEle)
        deltaHeight = 0;
        edgeDist = 0;
    elseif(deltaHeight<=0)
        deltaHeight = ventENU(3);
    end
    
%     ppENU = rENU
end

