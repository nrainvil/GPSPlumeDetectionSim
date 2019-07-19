function [ rv_az, rv_tang ] = findAzLim( vent_v, rcvr_v, rad)
%FINDAZLIM 
% Inputs:   vent_v: (ENU) Vent Location
%           rcvr_v: (ENU) Receiver location
%           rad: (m) Radius
% Outputs:  rv_az: (degrees) az limits 
%           rv_tang: (
rcvr_vent_v = vent_v - rcvr_v;

if(norm(rcvr_vent_v)<=rad)
    rv_az(1,:) = 0;
    rv_az(2,:) = 360;
    rv_tang(1,:) = [0,0];
    rv_tang(2,:) = [0,0];
    return;
end

rv_theta = acosd(rad/norm(rcvr_vent_v));

rv_tang(1,:) = rad.*(rcvr_vent_v(1:2)./norm(rcvr_vent_v(1:2)))*...
        [cosd(rv_theta),-sind(rv_theta);sind(rv_theta),cosd(rv_theta)];
rv_tang(2,:) = rad.*(rcvr_vent_v(1:2)./norm(rcvr_vent_v(1:2)))*...
        [cosd(-rv_theta),-sind(-rv_theta);sind(-rv_theta),cosd(-rv_theta)];
    
rv_tang(1,:) = vent_v(1:2) - rv_tang(1,:);    
rv_tang(2,:) = vent_v(1:2) - rv_tang(2,:);   

rcvr_tang_v(1,:) = rv_tang(1,:) - rcvr_v(1:2);
rcvr_tang_v(2,:) = rv_tang(2,:) - rcvr_v(1:2);
rv_az(1,:) = (180/pi())*atan2(rcvr_tang_v(1,2),rcvr_tang_v(1,1));
rv_az(2,:) = (180/pi())*atan2(rcvr_tang_v(2,2),rcvr_tang_v(2,1));
rv_az = rv_az+360;
end

