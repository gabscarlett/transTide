classdef RunConditions < handle
    % RUNCONDITIONS: Holds the run conditions for a tidal turbine
    % simulation. The geometry (blade profile) and environmental 
    % operating conditions are set
    
    % possible properties
    
    properties (Access = public)
               
        
        TurbineName (1,1) string = "none";
        
        Blades (1,1) {mustBeNumeric, mustBeInteger} = 3;
        
        RadialCoord (1,:){mustBeNumeric, mustBeFinite}= [] ; % radial coordinates
        
        BladeTwist (1,:){mustBeNumeric, mustBeFinite}= []; % geometrical blade twist along the radius (radians)
                        
        BladeChord (1,:){mustBeNumeric, mustBeFinite} = []; % chord length along the radius
        
        Depth (1,1){mustBeNumeric, mustBeFinite}= 45 ; % water depth
        
        HubDepth (1,1){mustBeNumeric, mustBeFinite}= -26 ; % depth of hub centre     
        
        TipSpeedRatio (1,1){mustBeNumeric, mustBeFinite}= 4.5 ; % ratio of flow to blade speed at the tip
        
        YawAngle (1,1){mustBeNumeric, mustBeFinite}= 0 ; % angle of rotor with onset flow
        
        PitchAngle (1,1){mustBeNumeric, mustBeFinite}= 0.017 ; % operational blade pitch angle (steady)
        
        HubVelocity (1,1){mustBeNumeric, mustBeFinite}= 2.7 ; % steady/mean flow velocity at the hub centre
        
        ShearLaw (1,1){mustBeNumeric, mustBeFinite}= 1/7 ; % exponent of shear velocity profile
        
        Waves = struct('On', 0, 'Height', [], 'Period', [], 'Direction', [], 'Type', [], 'Model', [], 'Periods', []);   % Structure containing the values for wave generation.
        
        Turbulence = struct('On', 0, 'Spectrum', [], 'Intensity', [], 'LengthScale', [], 'IsotropyRatio', []);   % Structure containing the values for synthetic turbulence.
        
        OperationsTab; % matlab table holding raw file data
        
        TurbTab;
        
       
    
    end
    
    methods
        
        function [obj] = RunConditions(varargin)
            
            % parse option to set run conditions from input file
            % for now we are assuming the file extension is .csv
            
            [opts, args] = checkOptions({{'turbine file',1},{'operating file',1}}, varargin);
            
            % check_options here
            
            if opts(1) % set tidal turbine parameters from file
                fileNameTurb = args{1}; % file name
                % open the file as a matlab table
                turbTab = readtable([fileNameTurb '.csv']);
                
                % set properties
                obj.TurbTab = turbTab;
                obj.TurbineName = turbTab.name{1};
                obj.Blades = turbTab.blades(1);
                obj.RadialCoord = turbTab.rad_coord;
                obj.BladeTwist = turbTab.twist;
                obj.BladeChord = turbTab.chord;
                
            end
            
            % operating conditions
            if opts(2) % set the operating conditions from file
                fileNameOps = args{2};
                % open the file as a matlab table
                opTab = readtable([fileNameOps '.csv']);
                
                % set properties
                obj.OperationsTab = opTab;                
                obj.Depth = opTab.depth(1);
                obj.HubDepth = opTab.hub_depth(1);
                obj.TipSpeedRatio = opTab.TSR(1);
                obj.YawAngle = opTab.yaw(1);
                obj.PitchAngle = opTab.pitch(1);
                obj.HubVelocity = opTab.hub_velocity(1);
                obj.Waves.On = opTab.waves(1);
                if obj.Waves.On
                    % set the struct
                    obj.Waves.Type = opTab.wave_type(1);
                    obj.Waves.Model = opTab.wave_model(1);
                    obj.Waves.Height = opTab.wave_H(1);
                    obj.Waves.Period = opTab.wave_T(1);
                    obj.Waves.Direction = opTab.wave_dir(1);
                    if strcmp(obj.Waves.Type, 'Irregular')
                        obj.Waves.Periods = opTab.wave_period_range(1):opTab.wave_period_range(2):opTab.wave_period_range(3);
                    end
                end
                obj.Turbulence.On = opTab.turbulence(1);
                if obj.Turbulence.On
                    % set the struct
                    obj.Turbulence.Spectrum = opTab.turb_spec{1};
                    obj.Turbulence.Intensity = opTab.turb_intensity(1);
                    obj.Turbulence.LengthScale = opTab.len_scale(1);
                    obj.Turbulence.IsotropyRatio = opTab.iso_ratio(1);      
                end
                
            end
 
        end
        
    end
    
end