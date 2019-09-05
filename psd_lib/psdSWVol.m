function [ sW_v, nmax ] = psdSWVol( DREVolume, k, uN, D_v )
%PSDSWVOL Return Particle Size Distributino based on 
% Dense Rock Equivalent Volume
% Input:    DREVolume
%           k: Shape factor
%           uN: Mode of particle size
%           D_v: Vector of diamter sizes
% Output:   sW_v: Vector of particle size counts

vol=DREVolume;

area_v = (4/3).*pi().*((D_v./2).^3);
volFcn = @(nm) sum(psdSW(k,uN,nm,D_v).*area_v);
volSearch = @(nm) volFcn(nm) - DREVolume;
nmax = fzero(volSearch,vol);

[ sW_v, ~, ~ ] = psdSW( k,uN,nmax,D_v);


end

