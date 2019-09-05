function [ lati,loni,alti ] = read_contour_map(whichVolcano)
% Input: whichVolcano - 1 for Etna, 2 for Redoubt
% Output: Array of lat,long,alt vectors representing DEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lati = 0;
loni = 0;
alti = 0;

if whichVolcano == 1
  %Load Data for 3-D Contour Map  
    fids = fopen('dem/Etna_Contour.csv','r');
elseif whichVolcano == 2
    fids = fopen('dem/redoubt_contour.csv','r');
else
    disp('Select Mt Etna (1) or Mt Redoubt (2)');
    return;
end

C = textscan(fids,'%f%f%f','Delimiter',',','Headerlines',1);
fclose(fids);
lati=[C{1}];
loni=[C{2}];
alti=[C{3}];


end

