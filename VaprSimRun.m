classdef VaprSimRun < handle
    %VAPRSIMRUN Batch run VaprSim at multiple times
    
    properties
        startTime = datetime(2018,1,1,1,0,0);  % (datetime) Current time 
        endTime = datetime(2018,1,1,1,1,0);    % (datetime) Current time
        deltaTime = duration(0,1,0);           % (duration) Simulation step
        sim;                                   % VaprSim object
        atn_c = {};                            % (dB) Simulation Results
        len_c = {};                            % (dB) Simulation Results
        prn_c = {};                            % (dB) Simulation Results
        frame_v = [];                          % animation frames
        paths_c = {};
        prns_c = {};
        cfg_plume = [];                        % Plume evolution matrix
        cfg_time = [];                         % Plume evolution matrix
    end
    
    methods
        function obj=VaprSimRun(startTime,endTime,vsim)            
            % VaprSimRun Constructor for VaprSimRun Object
            % Input:    startTime (datetime)
            %           endTime (datetime)
            %           vsim (VaprSim)
                obj.startTime = startTime;
            obj.endTime = endTime;
            obj.sim = vsim;
        end
        
        function runSim(obj,createGif)
            % runSim
            if(nargin<=1)
                createGif = 0;
            end
            
            dt_v = obj.startTime:obj.deltaTime:obj.endTime;
            obj.frame_v = [];
            
%             atn_lcl_c = {};
%             paths_lcl_c = {};
%             prns_lcl_c = {};
            for dtIdx = 1:length(dt_v)
                
                dt = dt_v(dtIdx);
                
                if(~isempty(obj.cfg_plume) && ~isempty(obj.cfg_time))
                    
                plumeLoc = copy(obj.sim.vents(1));
                if(dt<=obj.cfg_time(1))
                    cfg_plume_now = obj.cfg_plume(1,:);
                elseif(dt>=obj.cfg_time(end))
                    cfg_plume_now = obj.cfg_plume(end,:);
                else
                    time_comp = obj.cfg_time-dt;
                    [~,time_lb_loc] = min(abs(time_comp(time_comp<0)));
                    time_lb = obj.cfg_time(time_lb_loc);
                    time_ub_loc = time_lb_loc+1;
                    time_ub = obj.cfg_time(time_ub_loc);
                    time_scale = (dt-time_lb)/(time_ub-time_lb);
                    
                    cfg_plume_now = obj.cfg_plume(time_lb_loc,:) + ...
                             time_scale*(obj.cfg_plume(time_ub_loc,:) - ...
                                         obj.cfg_plume(time_lb_loc,:));
                end

                sA = cfg_plume_now(1,1);
                sB = cfg_plume_now(1,2);
                shiftDist = cfg_plume_now(1,3);
                shiftDir(1,1) = cfg_plume_now(1,4);
                shiftDir(1,2) = cfg_plume_now(1,5);
                plumeHeight = cfg_plume_now(1,6);
                MER = cfg_plume_now(1,7);
                Vvent= cfg_plume_now(1,8);
                
                shiftVec = shiftDist.*shiftDir./norm(shiftDir);
                plumeLoc.shiftENU([shiftVec,0],obj.sim.areaLoc);
                obj.sim.plumes(1).vent = plumeLoc;
                obj.sim.plumes(1).addMEREllCylPlume(...
                     sA,sB,shiftDir,plumeHeight,MER,Vvent,obj.sim.areaLoc);
                fprintf("Plume Size %0.0fx%0.0fm at %0.0fm\n",sA,sB,shiftDist);
                end 
                 
                [obj.atn_c{dtIdx}, obj.len_c{dtIdx}, ...
                    obj.prn_c{dtIdx}] = obj.simEpoch(dt,createGif);
                obj.paths_c{dtIdx} = obj.sim.signalPaths;
                obj.prns_c{dtIdx} = obj.sim.signalPRNs;
            end
%             obj.atn_c = atn_lcl_c;
%             obj.paths_c = paths_lcl_c;
%             obj.prns_c = prns_lcl_c;
        end
        
        function runMERSim(obj,mer,createGif)
            % runMERSim
            if(nargin<=1)
                createGif = 0;
            end
            
            dt_v = obj.startTime:obj.deltaTime:obj.endTime;
            plume = obj.sim.plumes(1);    
            plume.slices_v = {};
            obj.frame_v = [];
            
            for dtIdx = 1:length(dt_v)
                %update plume
                elapsedTime = dt_v(dtIdx)-obj.startTime;
                plume.updateMERSlices(mer,elapsedTime,...
                                      obj.deltaTime,obj.sim.areaLoc);
                
                dt = dt_v(dtIdx);
                [obj.atn_c{dtIdx}, obj.len_c{dtIdx}, ...
                    obj.prn_c{dtIdx}] = obj.simEpoch(dt,createGif);
            end
        end
        
        function [atn_rp_c, len_rp_c, prn_rp_c] = ...
                simEpoch(obj,dt,createGif)
           %simEpoch run one epoch
           ele = 10;
           fprintf('Running simulation at %s \n',dt);
           obj.sim.setTime(dt);
           obj.sim.findSignalPaths(ele);
           pFlag = 0;
           for rcvrIdx=1:length(obj.sim.receivers)
               atn_rp_v = [];
               len_rp_v = [];
               parfor pathIdx=1:length(obj.sim.signalPaths{rcvrIdx}(:))
                %for pathIdx=1:length(obj.sim.signalPaths{rcvrIdx}(:))
                   ray = obj.sim.signalPaths{rcvrIdx}(pathIdx);
                   prn = obj.sim.signalPRNs{rcvrIdx}(pathIdx);
                   [pathMass, pathLen, atn] = obj.sim.findSignalAtten(ray);
                   if(pathMass>0)
                       pFlag=1;
                       fprintf(' Intersection at %s (%d)',...
                           obj.sim.receivers(rcvrIdx).name,rcvrIdx);
                       fprintf(' PRN: %d (%d)',...
                           prn,pathIdx);
                       fprintf(' Length %0.2f (m)',pathLen);
                       fprintf(' Attenuation: %0.4f (dB)\n',atn);
                       % obj.sim.findSignalAtten(ray,1);
                   end
                   atn_rp_v(pathIdx) = atn;
                   len_rp_v(pathIdx) = pathLen;
                   prn_rp_v(pathIdx) = prn;
               end
               atn_rp_c{rcvrIdx} = atn_rp_v;
               len_rp_c{rcvrIdx} = len_rp_v;
               prn_rp_c{rcvrIdx} = prn_rp_v;
           end
           if(pFlag==0)
            fprintf('\n');
           end
           if(createGif==1)
               fig1 = obj.sim.plotSim(1);
               obj.frame_v = [obj.frame_v,getframe(fig1)];
               close(fig1);
           end
        end
        
        function fig_c = plotSim(obj,hideFig)
            
            if(nargin<=1)
                hideFig=0;
            end
            
            dt_v = obj.startTime:obj.deltaTime:obj.sim.time;
            for rcvrIdx=1:length(obj.sim.receivers)
                for pathIdx=1:length(obj.sim.signalPaths{rcvrIdx}(1,:))
                    for dtIdx=1:length(dt_v)
                        if(length(obj.atn_c{dtIdx}{rcvrIdx})>=pathIdx)
                            atn_v(dtIdx) = obj.atn_c{dtIdx}{rcvrIdx}(pathIdx);
                            len_v(dtIdx) = obj.len_c{dtIdx}{rcvrIdx}(pathIdx);
                        else
                            atn_v(dtIdx) = 0;
                            len_v(dtIdx) = 0;
                        end
                    end
                    
                    maxAtn = max(abs(atn_v));
                    if(max(atn_v)>0)
                        rcvr = obj.sim.receivers(rcvrIdx);
                        prn = obj.sim.signalPRNs{rcvrIdx}(pathIdx);
                        
                        if(hideFig==1)
                            fig = figure('Visible','off');
                        else
                            fig = figure;
                        end
                        
                        yyaxis left;
                        plot(dt_v,-1*atn_v,'-o','MarkerSize',4);
                        ylim([-maxAtn-.5,0]);
                        ylabel('Relative Signal Strength (dB)');
                        title(sprintf('%s PRN:%d',rcvr.name,prn));
                        yyaxis right;
                        plot(dt_v,len_v,'-o','MarkerSize',4);
                        fig_c{rcvrIdx}{pathIdx} = fig;
                    end
                end
            end
        end
        
        function isect_m = findIsect(obj)
            dt_v = obj.startTime:obj.deltaTime:obj.sim.time;
            for dtIdx=1:length(dt_v)
            for rcvrIdx=1:length(obj.sim.receivers)
                prnCnt = 1;
                for pathIdx=1:length(obj.sim.signalPaths{rcvrIdx}(1,:))
                    prn = obj.sim.signalPRNs{rcvrIdx}(pathIdx);
                    atn = obj.atn_c{dtIdx}{rcvrIdx}(pathIdx);
                    if(atn>0)
                        isect_m(dtIdx,rcvrIdx,prnCnt) = prn;
                        prnCnt = prnCnt+1;
                    end
                end
            end
            end
        end
    end
    
end

