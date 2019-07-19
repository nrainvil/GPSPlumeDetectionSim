function [ kappa_e ] = abs_coeff( freq, rho, e_s_rel, m)
%Absorption Coeffecient in dB/km
%
%rho - density (g/m^3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


u0 = 4e-7*pi();         %Vacuum permeability
e0 = 8.854187817e-12;   %Vacuum permittivity
e_s = e_s_rel.*e0;      %Scattering medium permittivity
e_s_pp = -1*imag(e_s);
w = 2*pi()*freq;
k = w.*sqrt(u0*e0);

kappa_e = k.*(m./rho).*(e_s_pp/e0).*abs(3*e0/(e_s+2*e0)).^2;

kappa_e = 4.34e3*kappa_e; %dB conversion

end

