function [sp3, found_orbits ] = load_sp3_information(year,month,day )
%function [sp3, found_orbits ] = load_sp3_information(year,month,day )
% inputs: year, month, day
% output: sp3 data and boolean (found_orbits)
found_orbits = false;
if year < 2000
  % going to assume you gave me a two character year
  year = year + 2000;
end
[GPS_wk, GPS_sec] = GPSweek(year,month,day,0,0,0);
gps_day = round(GPS_sec/86400);
% you could also make it go find the orbits, but wget does
% not currently work on my mac
% check for grm and then igs
% assume that the files live in a subdirectory called sp3
analysis_center='grm';
filename = ['sp3\' analysis_center num2str(GPS_wk) num2str(gps_day) '.sp3'];
analysis_center='igs';
filename2 = ['sp3\' analysis_center num2str(GPS_wk) num2str(gps_day) '.sp3'];
 
sp3 = read_sp3rev(filename);
if length(sp3) > 0
  found_orbits = true; 
else
  sp3 = read_sp3rev(filename2);
  if length(sp3) > 0
      found_orbits = true; 
  else
      fprintf("Missing %s\n", filename);
      fprintf("Missing %s\n", filename2);
    disp('no orbits, no simulation')
  end
end

end

