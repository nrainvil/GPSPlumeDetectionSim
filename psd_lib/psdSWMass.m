function [ sW_v, nmax ] = psdSWMass( DREMass, k, uN, D_v )
%PSDSWMASS Return Particle Size Distribution based on 
% Dense Rock Equivalent Mass
% Input:    DREMass (g)
%           k: Shape factor
%           uN: Mode of particle size
%           D_v: Vector of diamter sizes
% Output:   sW_v: Vector of particle size counts

rhoAsh = 2.813*(100^3); %g/m^3 Suwanosejima, Oguchi 2009 

area_v = (4/3).*pi().*((D_v./2).^3);
massFcn = @(nm) sum(psdSW(k,uN,nm,D_v).*area_v.*rhoAsh);
massSearch = @(nm) massFcn(nm) - DREMass;
nmax = fzero(massSearch,DREMass);

[ sW_v, ~, ~ ] = psdSW( k,uN,nmax,D_v);


end

