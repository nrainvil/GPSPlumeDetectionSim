function [ext_coeff] = getExtCoeff(obj,freqIn,cellX,cellY,cellZ)
%GETEXTCOEFF Return the extinction coeffecient for the selected cell
% Input:    cellX (pixel)
%           cellY (pixel)
%           cellZ (pixel)
% Output:   ext_coeff (N/m) Extinction Coeffecient

dens_c = obj.density(cellX,cellY,cellZ);
pdens_c = obj.particleDensity(cellX,cellY,cellZ);
dia_c = obj.diameter(cellX,cellY,cellZ);
prm_c = obj.permittivity(cellX,cellY,cellZ);
type_cnt = length(dens_c);

ext_coeff = 0;
for ii=1:type_cnt
    c_a = cross_abs(dia_c{ii}/2, freqIn, prm_c{ii});
    c_s = cross_sct(dia_c{ii}/2, freqIn, prm_c{ii});
    pmass = pdens_c{ii}.*((4/3).*pi().*(dia_c{ii}/2).^3);
    ndens = dens_c{ii}./pmass;
    ext_coeff = ext_coeff + ndens.*(c_a+c_s);
end

end

