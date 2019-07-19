function [ x_s ] = mie_cross_sct( a, freq, e_s_rel )
%MIE_CROSS_SCT 
% Mie Scattering Cross Section for a Spherical Particle
% Uses matscat library to find cross-section
%
% Input:    a - particle radius (m)
%           freq - frequency (Hz)
%           e_s_rel - Permittivity of Scattering medium (relative to e0)
% Output:   x_s - scattering cross section (m^2)
%
%%%%%%%%%%

%Run from VAPR_sim, VAPR_mie
addpath('../VAPR_mie/matscat');
addpath('../VAPR_mie/matscat/bessel');
addpath('../VAPR_mie/matscat/expcoeff');
addpath('../VAPR_mie/matscat/util');

x_s = zeros(size(a));
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
    x_s(kk) = C.sca;
end


end

