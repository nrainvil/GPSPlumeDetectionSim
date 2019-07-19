function [ obsCnt_v, cov, azOut, elOut ] = ...
    simSiteView( rcvrIdx, vent, sim, plumeRad, plumeHeight, dt_v )
%SIMSITEVIEW - Observations % per area
center = sim.areaLoc;
rcvr = sim.receivers(rcvrIdx);
azOut = zeros(length(dt_v),32);
elOut = zeros(length(dt_v),32);
for dtIdx=1:length(dt_v)
    dt = dt_v(dtIdx);
    sim.setTime(dt);
    txt_c{dtIdx} = {};
    
    [rENU(1), rENU(2), rENU(3)] = rcvr.enu(center);
    %Plot Vents
    [vENU(1), vENU(2), vENU(3)] = vent.enu(center);

    for pIdx=1:length(sim.signalPaths{rcvrIdx})
        path = sim.signalPaths{rcvrIdx}(pIdx);
        prn = sim.signalPRNs{rcvrIdx}(pIdx);
        [azel(1),azel(2)] = path.getAzEl();
        azel(1) = gpsAz(azel(1));
        azOut(dtIdx,prn) = azel(1);
        elOut(dtIdx,prn) = azel(2);
        %                 fprintf('%s PRN:%d Az: %0.2f El: %0.2f\n',...
        %                     rcvr.name,prn,azel(1),azel(2));
        [edge_dist, delta_height, ppENU ] = ...
            findPP( vent, rcvr, azel, plumeRad, plumeHeight, center);
        if(edge_dist>0)
            txt_c{dtIdx}{end+1} = {rcvr.name,prn};
            az_adj = matAz(azel(1));
            az_ref_v = edge_dist*[cosd(az_adj),sind(az_adj)]; %Direction vector
%             fprintf('%s intersect on PRN:%d',rcvr.name,prn);
%             fprintf(' Dist: %0.2f m Height: %0.2f m El: %0.2f\n',...
%                 edge_dist,delta_height,elOut(dtIdx,prn));
        else
            elOut(dtIdx,prn) = 0;
        end
    end

end

%%
cov=0;
obsCnt_v = zeros(size(dt_v));
for dtIdx=1:length(dt_v)
    dt = dt_v(dtIdx);
    if(~isempty(txt_c{dtIdx})>0)
%         fprintf('%02.0f:%02.0f',dt.Hour,dt.Minute);
        txtLen = length(txt_c{dtIdx});
        obsCnt_v(dtIdx) = txtLen;
%         for txtIdx=1:txtLen
%             fprintf(' %s:%d',txt_c{dtIdx}{txtIdx}{1},txt_c{dtIdx}{txtIdx}{2});
%         end
        cov = cov + 1;
%         fprintf('\n');
    end
end
cov = 100.*cov./length(dt_v);
end

