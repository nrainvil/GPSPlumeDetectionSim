function [sp3] = read_sp3rev( filename )
%------------------------------------------------------------------------%
% This code reads Cartesian satellite positions from an SP3 file
%
% Inputs:   filename - filename of SP3 file
%
% Outputs:  sp3_obs = [GPS_week GPS_TOW PRN x y z], sorted by PRN #
%           *note - GPS_TOW is in seconds;  x,y,z are in km
%	     clockbias unit is microseconds
% sp3_obs is a structured variable, data and col.
% col defines what is in each column. But you  
% will see that gps week and seconds are in columns 1 and 2.  PRN is 
% in column 3, xyz in columns 4-6, and clockbias in 7.
% REV:
% column 8 now has the GNSS constellation type

%
% Coded by David Wiese
% Colorado Center for Astrodynamics Research 
% University of Colorado at Boulder
% October 12, 2006

% 15sep18
% modified by Kristine Larson to return the GNSS constellation 
% type in column 8
% 1 GPS
% 2 GLONASS
% 3 GALILEO
% 4 COMPASS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp3 = [];
 
% Assign a file ID and open the given header file.
fid=fopen(filename);

% If the file does not exist, display warning message
if fid == -1
%   display(['SP3 file does not exist: ' filename]);
  return
else
%   display(['Using: ' filename]);
     
  % Go through the header (23 lines long)
  for i = 1:23
    current_line = fgetl(fid);
         
    % Store the number of satellites in the SP3 file    
    if i == 3
      current_line = current_line(2:length(current_line));
      F = sscanf(current_line,'%u');
      no_sat = F(1);
    end
    i = i + 1;
  end
    
    % Begin going through times and observations
  end_of_file = 0;
  i = 0; j = 1;
  while end_of_file ~= 1        
    current_line = current_line(2:length(current_line)) ;
    F = sscanf(current_line,'%f');
    % Load GPS Gregorian time into variables
    Y = F(1);
    M = F(2);
    D = F(3);
    H = F(4);
    min = F(5);
    sec = F(6);
    Greg_time(j,:) = [Y M D H min sec];
        
     % Convert GPS Gregorian time to GPS week and GPS TOW
    [GPS_wk, GPS_TOW] = GPSweek(Y,M,D,H,min,sec) ;
         
     % Store satellite PRN and appropriate observations
    for n = 1:no_sat
       % Go to the next line
       current_line = fgetl(fid);
       % figure out the constellation type
       constellation =  constellation_type(current_line(2:2));
       current_line = current_line(3:length(current_line));
       F = sscanf(current_line,'%f');
            
            % Save PRN, positions, and clock error
       PRN = F(1); x = F(2); y = F(3); z = F(4); clk_err = F(5);
            
            % Create observation vector
            % KL augment with constellation type (1-4)
       sp3_obs_all(i+n,:) = [GPS_wk GPS_TOW PRN x y z clk_err constellation];
       n = n + 1;
    end       
        % Go to next line - check to see if it is the end of file
    current_line = fgetl(fid);
    if strfind(current_line,'EOF')
      end_of_file = 1;
    end
        
    i = i + n - 1;
    j = j + 1;
  end         
end

sp3.data = sp3_obs_all;
sp3.col.WEEK = 1;
sp3.col.TOW = 2;
sp3.col.PRN = 3;
sp3.col.X = 4;
sp3.col.Y = 5;
sp3.col.Z = 6;
sp3.col.B = 7;
sp3.col.GNSS = 8;



