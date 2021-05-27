function [Cd_DS_3D] = dynStallCd(Dvis, Cd_Ind, grid, aE, r, Cd0)

if size(grid.FCd.GridVectors,2) == 2
    R_InRange=(r>0.8*r(end));        % boolean (TRUE for outer section)
    % DRAG 3D
    
    % complete unsteady drag coefficient
    
    % No rotational correction for outer sections to avoid over
    % prediction (drag reduces from 0.8R - R during rotation)
    Cd_2St = grid.Cd2d(aE(R_InRange,:));
    rr=r(~R_InRange)'.*ones(size(aE(~R_InRange,:)));
    
    %%
    
    Cd_3St = grid.FCd(aE(~R_InRange,:), rr);
    Cd_St = vertcat(Cd_3St, Cd_2St);
    
else
    
    Cd_St = grid.FCd(aE);
end

% complete unsteady drag coefficient
Cd_DS_3D = Cd_Ind + Cd_St + (Cd_St-Cd0).*Dvis;

end