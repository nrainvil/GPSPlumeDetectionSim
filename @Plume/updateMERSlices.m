function updateMERSlices(obj, MER, elapsedTime, deltaTime, centerLocation)
%UPDATEMERPLUME Createa a plume based on the MER
% Inputs:       mer  (kg/s) Mass eruption rate
%               elapsedTime (duration)

gAccel = 9.8; % (m/s^2)
ventVel = [0,0,184]; % (m/s) Feret-Lorgeril
ventDia = 100; %(m)
[ventEnu(1),ventEnu(2),ventEnu(3)] = obj.vent.enu(centerLocation); %(m)
tsD = duration(0,0,1); % s
ts = seconds(tsD);
dtMass = MER*seconds(ts); %(kg)

startTime = max([elapsedTime-deltaTime,duration(0,0,0)]);    
dt_v = startTime:tsD:elapsedTime;
fprintf('Update MER: %s\n',elapsedTime);

for dtIdx = 1:length(dt_v)
    dt = seconds(dt_v(dtIdx)); % s
    fprintf('%d|%d ',dtIdx,dt);
    
    for slRevIdx = 1:length(obj.slices_v)+1
        %Assumes acceleration in vertical direction only
        %Assume constant height slice
        slIdx = length(obj.slices_v) - slRevIdx + 2;
        fprintf('%d,',slIdx);
        if(slIdx==1)
            sliceHeight = ventVel(3)*ts; % (m)
            sliceVel = [ventVel(1),ventVel(2),ventVel(3) - gAccel*ts];% m/s
            sliceENU = [ventEnu(1) + ventVel(1)*ts, ...
                ventEnu(2) + ventVel(2)*ts, ...
                ventEnu(3) + ((ventVel(3)+sliceVel(3))/2)*ts];
            obj.slices_v{slIdx} = ...
                sliceMER(ventDia,sliceHeight,sliceENU,dtMass,sliceVel);
        else
            obj.slices_v{slIdx} = obj.slices_v{slIdx-1};
            oldENU =  obj.slices_v{slIdx}.enu;
            oldVel =  obj.slices_v{slIdx}.velocity;
            newVel = [oldVel(1),oldVel(2),oldVel(3) - gAccel*ts];% m/s
            obj.slices_v{slIdx}.velocity = newVel;
            obj.slices_v{slIdx}.enu = ...
                [oldENU(1) + oldVel(1)*ts, ...
                oldENU(2) + oldVel(2)*ts, ...
                oldENU(3) + ((oldVel(3)+newVel(3))/2)*ts];
        end
    end
    fprintf('\n');  
end

plumeSize = size(obj.density);
fprintf('Updating Plume %dx%dx%d %d\n',plumeSize(1),plumeSize(2),...
                                       plumeSize(3),plumeSize(4));

%Iterate over density matrix
plumePS = obj.getPixSize();
for alt_idx=1:plumePS(3)
    for x_idx=1:plumePS(1)
        for y_idx=1:plumePS(2)
            
            obj.density(x_idx,y_idx,alt_idx,:) = 0;
        end
    end
end

end

