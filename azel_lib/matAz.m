function [ azMat ] = matAz( azIn )
%MATAZ 
while(azIn<0)
    azIn=mod(azIn+360,360);
end

azMat = (360-azIn)+90;
azMat = mod(azMat,360);

end

