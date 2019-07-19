function [ x_s ] = cross_sct( a, freq, e_s_rel )
%CROSS_SCT 
% Raleigh Scattering Cross Section for a Spherical Particle
%
% Input:    a - particle radius (m)
%           freq - frequency (Hz)
%           e_s_rel - Permittivity of Scattering medium (relative to e0)
% Output:   x_s - scattering cross section (m^2)
%
%%%%%%%%%%

u0 = 4e-7*pi();         %(H/m) Vacuum permeability
e0 = 8.854187817e-12;   %(F/m) Vacuum permittivity
e_s = e_s_rel*e0;       %Scattering permittivity

w = 2*pi()*freq;        % (rad/s)
k = w.*sqrt(u0*e0);     % (rad/m) -- rad/s*sqrt(s^2/m^2)

x_s = (8*pi()/3).*k.^4.*a.^6.*abs((e_s-e0)./(e_s+2*e0)).^2;

end

