function [log_m] = bin_deg(deg_m,deg,deg_bin)
%BIN_DEG
% Return logical array containing values where deg_m is near deg

deg_min = mod(deg - deg_bin/2,360);
deg_max = mod(deg + deg_bin/2,360);

if(deg_max>deg_min)
    log_m = (deg_m > deg_min) & (deg_m < deg_max);
else
    log_m = (deg_m < deg_max) | (deg_m > deg_min);
end

end

