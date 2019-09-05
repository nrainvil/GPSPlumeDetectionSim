function [ azGPS ] = gpsAz( azIn )
%GPSAZ 
if(diff(azIn)==360)
    azGPS = [0;360];
    return;
end

while(azIn<0)
    azIn=mod(azIn+360,360);
end

azGPS = (360-azIn)+90;
azGPS = mod(azGPS,360);

end

