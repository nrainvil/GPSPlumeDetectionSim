function [ out ] = constellation_type( onechar )
%function [ out ] = constellation_type( onechar )
% send it one character (G,R,E, etc) from a sp3 file (or RINEX)
% and it returns an integer representation for the constellation
% 1 GPS
% 2 GLONASS
% 3 GALILEO
% 4 COMPASS
% Author: Kristine Larson
% 15sep17

% default is GPS
out = 1;

switch onechar
    case 'G'
        out=1;
    case 'R'
        out=2;
    case 'E';
        out=3;
    case 'B';
        out=4;  
end

end

