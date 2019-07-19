function [ x_a ] = cross_abs( a, freq, e_s_rel )
%CROSS_ABS 
% Raleigh Absorption Cross Section for a Spherical Particle
%
% Input:    a - particle radius (m)
%           freq - frequency (Hz)
%           e_s_rel - Permittivity of Scattering medium (relative to e0)
% Output:   x_a - absorption cross section (m)
%
%%%%%%%%%%

u0 = 4e-7*pi();         %Vacuum permeability
e0 = 8.854187817e-12;   %Vacuum permittivity
e_s = e_s_rel*e0;       %Scattering permittivity

e_s_pp = -1*imag(e_s);
w = 2*pi()*freq;
k = w.*sqrt(u0*e0);     %(rad/m) -- rad/s*sqrt(s^2/m^2)

x_a = (4*pi()/3)*k.*a.^3.*(e_s_pp/e0).*abs(3*e0/(e_s+2*e0)).^2;

end

