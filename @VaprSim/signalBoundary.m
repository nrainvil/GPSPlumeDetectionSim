function [ coordsB ] = signalBoundary(obj, rayIn)
%SIGNALBOUNDARY Calculate boundry of RF signal relative to simulation 
%               environment
% Input:    signalIn (RFSignal)
% Output:   e (m) East
%           n (m) North
%           u (m) Up

pixSize_v = obj.getPixSize();

bPlane_v(1) = Plane3D([0,0,pixSize_v(3)],[0,0,1]); %up
bPlane_v(2) = Plane3D([0,0,0],[0,0,-1]); % down
bPlane_v(3) = Plane3D([pixSize_v(1)/2,0,0],[1,0,0]); % east
bPlane_v(4) = Plane3D([-1*pixSize_v(1)/2,0,0],[-1,0,0]); % west
bPlane_v(5) = Plane3D([0,pixSize_v(2)/2,0],[0,1,0]); % North
bPlane_v(6) = Plane3D([0,-1*pixSize_v(2)/2,0],[0,-1,0]); % South

for ii=1:length(bPlane_v)
    dist_v(ii) = bPlane_v(ii).rayIntersect(rayIn); 
end

[minDist,minLoc] = min(dist_v);
coordsB = rayIn.getXYZ(minDist).*obj.areaRes;
% coordsB = (minDist.*rayIn.dir + rayIn.origin).*obj.areaRes;

end

