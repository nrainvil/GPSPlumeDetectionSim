function [az,el] = ecef2azelrange(r_sat,r_site,latgd,lon)

%==========================================================================
%==========================================================================
% [az,el,range] = ecef2azelrange(r_sat,r_site,latgd,lon)
%
% Calculates the azimuth, elevation, and range of a satellite with respect
%  to an observation site.
%
%
% Author: Ben K. Bradley
% Date: 11/15/2010
% Modified to remove calculations for ASEN5090 assignments
%
% INPUT:         Description                                         Units
%
%  r_sat      - position of satellite in ECEF frame                 [x y z]
%  r_site     - position of observing site in ECEF frame            [x y z]
%  latgd      - geodetic latitude of observation site          [-90,90] deg
%  lon        - longitude of observation site     [-180,180] or [0,360] deg
%
%
% OUTPUT:       
%    
%  az         - azimuth (degrees clockwise from North)          [0,360] deg
%  el         - elevation (degrees up from horizon)            [-90,90] deg
 %
%
% Coupling:
%
%  none
%
%
%==========================================================================
%==========================================================================

% Processig Input
latgd_rad = degtorad(latgd);
lon_rad = degtorad(lon);
rho_ecef = r_sat - r_site;
rho_sez = rot2mat(pi/2-latgd_rad)* rot3mat(lon_rad) * rho_ecef;
rho = norm(rho_sez);

% Output
el = asind(rho_sez(3)/rho);
az = atan2d(rho_sez(2),-rho_sez(1));
if (az<0)
    az = az+360;
end
range = rho;

