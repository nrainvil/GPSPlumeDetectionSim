function [ x_a ] = mie_cross_abs( a, freq, e_s_rel )
%MIE_CROSS_ABS 
% Mie Absorption Cross Section for a Spherical Particle
%
% Input:    a - particle radius (m)
%           freq - frequency (Hz)
%           e_s_rel - Permittivity of Scattering medium (relative to e0)
% Output:   x_a - absorption cross section (m)
%
%%%%%%%%%%

%Run from VAPR_sim, MatScat
addpath('../MatScat');
addpath('../MatScat/bessel');
addpath('../MatScat/expcoeff');
addpath('../MatScat/util');

x_a = zeros(size(a));
lambda = physconst('LightSpeed')/freq; 

%Calculate Mie cross-sections
nang = 1800;        % number of far field angles to evaluate
conv = 1;           % convergence factor
nm = 1;             % outer medium refractive index (real)
ns = conj(sqrt(e_s_rel));   % ??? conjugate sphere refractive index (complex)
for kk=1:length(a)
    dia = 2*a(kk); %m
    rad = dia/2;           % sphere radius
    [~, C, ~] = calcmie(rad, ns, nm, lambda, nang, ...
        'ConvergenceFactor', conv);
    x_a(kk) = C.abs;
end

end

