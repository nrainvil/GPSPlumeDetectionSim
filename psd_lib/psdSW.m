function [ s_w, f_w, gamma ] = psdSW( k,u_n,nmax,D_v)
%PSDSW Scaled Weibull PSD
% Input:    k = shape factor
%           u_n = mode of PSD
%           nmax = scale factor
%           D_v = diameter vector
% Output:   s_w = PSD
%           f_w = PDF

gamma = u_n.*((k-1)/k).^(-1/k);

f_w = (k/gamma) .* ((D_v./gamma).^(k-1)) .* exp(-1*((D_v/gamma).^k));

s_w = (f_w/max(f_w)).*nmax;


% Dn = 7.25*u_n;
% u = k-1;
% v = k;
% Gn = gamma;
% 
% s_w = ((D_v./Dn).^u).*exp(-Gn.*(D_v./Dn).^v);
% s_w = nmax*s_w./max(s_w);


end

