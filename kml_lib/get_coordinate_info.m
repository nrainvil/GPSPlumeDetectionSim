function [x,y, x_stn, y_stn, XYZr] = get_coordinate_info(lati,loni,...
    alti,lat,lon,alt, station,volcano)
% inputs are 
%     the lat,long, and altitude of the DEM in the variables
%      called lati,loni,alti.
%      the lat, long, and altitude of the stations in the variables
%      called lat,lon,alt
% outputs
% x and y are the planar coordinates of the DEM, in km.
% x_stn and y_stn are the planar coordinates of the stations, in km
% these are with respect to the volcano
% XYZr are cartesian coordinates for the stations N by 3, where N is number
% of stations. units are meters
   %% 1 degree of latitude in meters at equator.
    Eq_lat = 111325; 
  % 1 degree of latitiude in meters  
    onedeglat = Eq_lat;
    onedeglon = cosd(volcano.lat) * Eq_lat;
  % Convert to local planar system.
    % x - reference
    volcano_y = volcano.lat * onedeglat; % x = 0;
    % y - reference  
    volcano_x = volcano.lon * onedeglon; % y = 0;
  % Convert Contour map of Lat Lon into (horizontal) planar coordinates 
    for i = 1:length(lati)     
      % Units (km) w.r.t volcano location
      y(i) = (lati(i)*onedeglat - volcano_y)/1000;
      x(i) = (loni(i)*onedeglon - volcano_x)/1000;
    end
  % Convert Station lat lon to (horizontal) planar coordinates 
    for i = 1:length(lat)      
      y_stn(i) = (lat(i)*onedeglat - volcano_y)/1000;
      x_stn(i) = (lon(i)*onedeglon - volcano_x)/1000;
      dd= sqrt(y_stn(i)*y_stn(i) + x_stn(i)*x_stn(i));
      fprintf('X %6.2f Y %6.2f %s (km) \n', x_stn(i), y_stn(i), station(i,1:4));
    % text(x_stn(i),y_stn(i),alt(i)+5,station(i,:),'FontSize',16);
      % real cartesian
    end
    %lat and lon are in degrees, altitude is in meters
    % now compute Cartesian receiver coordinates in meters?
    for i = 1:length(lat)      
      [xx,yy,zz]=lla2ecef(pi*lat(i)/180,pi*lon(i)/180,alt(i));
      XYZr(i,1:3) = [xx yy zz];
    end  


end

