function [ x_s] = cross_sct_jvgr( a, freq, e_s_rel )
%CROSS_SCT_JVGR 
% Raleigh Absorption Cross Section for a Spherical Particle
% using JVGR definition
%
% Input:    a - particle radius (m)
%           freq - frequency (Hz)
%           e_s_rel - Permittivity of Scattering medium (relative to e0)
% Output:   x_s - scattering cross section (m)
%
%%%%%%%%%%

u0 = 4e-7*pi();         %Vacuum permeability
e0 = 8.854187817e-12;   %Vacuum permittivity

w = 2*pi()*freq;
k = w.*sqrt(u0*e0);
c = 1/sqrt(u0*e0);
lambda = c/freq;
D = 2*a;

chi = pi().*D./lambda;
K = (e_s_rel-1)./(e_s_rel+2);
K2 = abs(K).^2;
x_s = (2*lambda.^2)/(3*pi()).*chi.^6.*K2;

end

