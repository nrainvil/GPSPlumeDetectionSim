function [ sat,XYZs,satid ] = do_orbits(sp3, year,month,day,Hour,Minute,maxConst )
% inputs
% year, month, day, hour, minute
% maxConst is the maximum number of constellations
% outputs
% vector array of satellites, their cartesian coordinates (meters?), 
% and satid (1-GPS,2-GLONASS,3-GALILEO)
% fprintf('%d:%d:%d %d:%d:%d\n',year,month,day,Hour,Minute,0);
maxsat = 32 ;% (only 32 satellites per constellation
[GPS_wk, GPS_sec] = GPSweek(year,month,day,Hour,Minute,0);
sat = []; XYZs = [];satid=[]; % arrays to store satellite names and coordinates at epoch
for constellation = 1:maxConst
  for s=1:maxsat       
    [satPos] = preciseO(GPS_sec, s, sp3,constellation) ;
    if length(satPos) > 0
     % save satellite name and coordinates (in meters)
       sat = [sat; s];
       XYZs = [XYZs; 1000*satPos];
       satid = [satid; constellation];
    end
  end % satellite loop
end % constellation loop

