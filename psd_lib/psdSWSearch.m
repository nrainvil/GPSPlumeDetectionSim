function [ vol ] = psdSWSearch( k,u_n,nmax,D_v)
%PSDSWSEARCH 

psdSWfun = @(D) psdSW(k,u_n,nmax,D);
volTot = sum(s_w.*(4/3).*pi().*((D_v./2).^3));

end

