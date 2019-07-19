function [ in ] = inEllipse( point,saA,saB,dirA )
%INELLIPSE Checks if 2D point is inside an ellipse
%          Assumes ellipse is centered at 0
% Input:    point (2x1)
%           saA (scalar) semi-axis A
%           saB (scalar) semi-axis B
%           dirA (2x1) unit vector in direction of A


dirA = dirA./norm(dirA);
RA = [dirA(1),dirA(2);-dirA(2),dirA(1)];
pointR = RA*point';
testRes = (pointR(1,:).^2)./(saA.^2) + (pointR(2,:).^2)./(saB.^2);
in = (testRes <= 1);

end

