function [ x_s ] = chrg_cross_sct( a, freq, e_s_rel, eta )
%MIE_CROSS_SCT 
% Charged Surface Mie Scattering Cross Section for a Spherical Particle
% Load look up table from csv
% Table calculated with modified PyMieScatt 
%
% Input:    a - particle radius (m)
%           freq - frequency (Hz)
%           e_s_rel - Permittivity of Scattering medium (relative to e0)
% Output:   x_s - scattering cross section (m^2)
%
%%%%%%%%%%
global g_cm_csv
if(isempty(g_cm_csv))
    g_cm_csv = csvread('vsr_chrg.csv');
end

x_s = [];
for kk=1:length(a)
    dia = 2*a(kk);
    cm_new = g_cm_csv(round(g_cm_csv(:,1),5)==round(dia,5) & ...
                g_cm_csv(:,2)==freq & ...
                g_cm_csv(:,3)==real(e_s_rel) & ...
                g_cm_csv(:,4)==abs(imag(e_s_rel)) & ...
                round(g_cm_csv(:,10),2,'significant')==round(eta,2,'significant'),:);
    if(isempty(cm_new))
        fprintf('Missing CS: %0.8f, %0.8f, %0.8f, %0.8f\n',dia,freq,real(e_s_rel),imag(e_s_rel));
        return;
    end
    x_s(kk) = cm_new(1,9)*(1e-9^2); %nm^2 to m^2
end


end

