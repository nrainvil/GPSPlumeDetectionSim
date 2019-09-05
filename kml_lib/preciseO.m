function [XYZpos] = preciseO(Ttime, PRN, sp3,constel)
%function [XYZpos] = preciseO(Ttime, PRN, sp3)
% returns Cartesian position of satellites (in meters)
% input is 
%     Time in GPS seconds, that you want the orbit
%     PRN number,
%     pre-read values for sp3 file
%    
%
% week is 1st, GPS seconds in 2nd
% PRN is in 3rd colum
% XYZ in 4-6
% you want the orbit location at Ttime
% 
% constants

% default value
XYZpos = [];
c=299792458;    % speed of light m/s
debug = true;
relativity = 0;
sidereal_day = 0.99726956634;
period = sidereal_day;
P0 = 2.0*pi / period;
T0 = Ttime;
ii = find(sp3.data(:,3) == PRN & sp3.data(:,8) == constel);
if length(ii) > 0
  OrigTime = sp3.data(ii,2);
  X = sp3.data(ii,4);
  Y = sp3.data(ii,5);
  Z = sp3.data(ii,6);
%find the nearest zero point
  [xx, closestI] = min(abs(Ttime- OrigTime));
  NDAT = 4;
% indices to use in the fit
if(closestI-NDAT>0)
  k = [closestI-NDAT:closestI+NDAT];
else
  k = 1:9;
end

%orig = find(OrigTime(k) == T0)
  Time = (sp3.data(ii,2)-T0)/86400;
  if(max(k)>length(Time))
      k=length(Time)-8:length(Time);
  end
%Time(k)

  Nest = 9;
  Nmeas = length(k);
  A  = zeros(Nmeas, Nest);
  A(:,1) = ones(Nmeas,1);  % needs to be sized
  B  = zeros(1, Nest);
  B(1,1) = 1;

gps_rel_time = 0;

gps_rel_time = (Ttime -T0)/86400;
ND = (Nest-1)/2;
for ii = 1:ND
  kk = 2+(ii-1)*2;
  P =  P0*ii;

  A(:,kk) = sin(P*Time(k));
  A(:,kk+1) = cos(P*Time(k));
  % this is stupid - define zero and one instead of using cosine/sine
  B(1,kk) = sin(P*gps_rel_time);
  B(1,kk+1) = cos(P*gps_rel_time);
end
%if debug
  % 86400*Time(k) + T0
%end
XCoeffs = A\X(k); newX = A*XCoeffs;
YCoeffs = A\Y(k); newY = A*YCoeffs;
ZCoeffs = A\Z(k); newZ = A*ZCoeffs;
%these are in meters
size(newX);
XYZpos = [B*XCoeffs B*YCoeffs B*ZCoeffs];
 
end

