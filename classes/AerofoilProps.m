classdef AerofoilProps < handle
    % AEROFOILPROPS: Holds the static aerofoil section properties for a simulation.
    % Given the coefficients for lift and drag with angle of attack 
    % the class computes normal and chordwise coefficients or vice verca 
    % computes the static seperation point and extrapolates values to 360 deg. 

    properties (Access = public)
               
        
        AerofoilName (1,1) string = "none";
        
        FoilTabData; % matlab table holding raw file data.
         
        ZeroLiftAngle (1,1){mustBeNumeric, mustBeFinite} = 0; % angle of attack corresponding to Cl = 0;
        
        ZeroLiftDrag (1,1){mustBeNumeric, mustBeFinite} = 0; % drag coefficient for alpha(Cl = 0)

        LinearClRange (1,:){mustBeNumeric, mustBeFinite} = deg2rad([-5, 5]) % range where Cl is linear with angle
        
        LinearLiftSlope (1,1){mustBeNumeric, mustBeFinite}= 2*pi; % slope in the linear Cl range
        
        AngleOfAttack (1,:){mustBeNumeric, mustBeFinite}= deg2rad(-5:0.1:5); % always in radians !!!
                        
        LiftCoeff (1,:){mustBeNumeric, mustBeFinite} = []; % lift coefficient
        
        DragCoeff (1,:){mustBeNumeric, mustBeFinite} = []; % drag coefficient
        
        NormCoeff (1,:){mustBeNumeric, mustBeFinite} = []; % normal coefficient
        
        AxialCoeff (1,:){mustBeNumeric, mustBeFinite}= []; % chordwise coefficient
    
    end
    
    properties (SetAccess = private)
        
        CoeffValues_360 = [];
    end
    
    methods
        
        function [obj] = AerofoilProps(fileName)
            
            foilTab = readtable([fileName '.csv']);
            
            obj.FoilTabData = foilTab;
            obj.AerofoilName = foilTab.name{1};
            obj.ZeroLiftAngle = foilTab.aoa_Cl0(1);
            obj.ZeroLiftDrag = foilTab.Cd_Cl0(1);
            obj.LinearClRange = [foilTab.Cl_lin_range(1) foilTab.Cl_lin_range(2)];
            obj.LinearLiftSlope = foilTab.Cl_lin_slope(1);
            
            obj.AngleOfAttack = foilTab.aoa;
            obj.LiftCoeff = foilTab.Cl;
            obj.DragCoeff = foilTab.Cd;
            
            % normal coefficient
            if ~isnan(foilTab.Cn(1))
                obj.NormCoeff = foilTab.Cn;               
            else
                obj.NormCoeff = obj.LiftCoeff.*cos(obj.AngleOfAttack)...
                    + obj.DragCoeff.*sin(obj.AngleOfAttack); % normal coefficient
            end
            
            % set coefficients for the full 360 degree range
            obj.setDeepStallValues;
            % set the static seperation point
            obj.setSeperationPoint;
 
        end
        
        function [] = setDeepStallValues(obj)
            
        % Deep stall extrapolation (-180 <-> 180)
        [Values_360.Alpha, Values_360.Cl, Values_360.Cd, Values_360.Cn] = vitExtrapolation(obj.AngleOfAttack, obj.LiftCoeff, obj.DragCoeff);
        obj.CoeffValues_360 = Values_360;
        
            
        end
        
        function [] = setSeperationPoint(obj)
            
            % Non-rotational separation point (-180 <-> 180)
            obj.CoeffValues_360.F = sepPoint(obj.CoeffValues_360.Alpha, obj.ZeroLiftAngle, obj.CoeffValues_360.Cn, ...
                obj.LinearLiftSlope,obj.LinearClRange); 
            
        end
        
    end
    
end