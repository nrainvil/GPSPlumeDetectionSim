function [ x_a ] = cross_abs_jvgr( a, freq, e_s_rel )
%CROSS_ABS_JVGR 
% Raleigh Absorption Cross Section for a Spherical Particle
% using JVGR definition
%
% Input:    a - particle radius (m)
%           freq - frequency (Hz)
%           e_s_rel - Permittivity of Scattering medium (relative to e0)
% Output:   x_a - absorption cross section (m)
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
x_a = (lambda.^2./pi()).*chi.^3.*imag(K);
x_a = abs(x_a); %??

end

