classdef TidalSim < handle
    % TIDALSIM: Runs a tidal turbine simulation for a given
    % geometry and operating conditions.
    % Computes flow, loads, bending moments and power.
    
    % possible properties
    
    properties (Access = public)
        
        BladeSections (1,1) {mustBeNumeric, mustBeInteger} = 100; % discretisation along the blade span
        Rotations (1,1) {mustBeNumeric, mustBeInteger} = 100; % number of rotor revolutions
        Steps (1,1){mustBeNumeric, mustBeInteger} = 72; % number of time steps per revolution
        Run; % class containing run conditions for the simulation
        SeedTurbulence = []; % random seed for turbulence model
        SeedWaves = []; % random seed for wave model
        AeroFoil; % class containing aerofoil static coefficients for the simulation
        RotationalAugmentation(1,1){mustBeNonnegative, mustBeNumericOrLogical, mustBeInteger, mustBeLessThanOrEqual(RotationalAugmentation,1)} = false;
        LoadMethod (1,1) string {mustBeMember(LoadMethod, ["Quasi-steady","Unsteady"])} = 'Quasi-steady'; % the loads can be computed by quasi-steady or unsteady methods
        Density = 1025; % fluid density
        DSData; % file path to empirical data for dynamic stall model (LoadMethod = 'Unsteady')
        
        
    end
    
    properties(SetAccess = private)
        Omega % rotational speed (rad/s)
        RotationPeriod % period of rotation (s)
        Phase % phase difference from leading blade
        Psi % temporal position of leading blade
        Time; % simulation time array
        RadialCoords; % radial coordinates after blade section discretisation
        BladeTwist; % geometrical twist and blade pitch after blade section discretisation
        BladeChord; % chord lengths after blade section discretisation
        Radius; % radius of rotor
        UAxial; % Axial velocity seen by the blade
        UTangential; % Tangential velocity seen by the blade
        URelative; % velative velocity at the blade
        AngleOfAttack; % Angle of attack (angle between velocity components)
        FlowAngle; % flow angle (angle of attack including blade twist/pitch angle)
        AeroCoeffs; %
        AxialIndFactor;  % axial induction factor computed using BEMT
        TangentialIndFactor % tangential induction factor computed using BEMT
        LiftCoeff; % Lift coefficient time series
        DragCoeff; % Drag coefficient time series
        ForceNormal; % Normal force at each blade section (N/m)
        ForceTangential; % Tangential force at each blade section (N/m)
        RootBM; % Root bending moment time series for each blade
        EdgeBM; % Edgewise bending moment time series for each blade
        Power; % Power time series for each blade
        Thrust; % Thrust time series for each blade
    end
    
    properties(Access = private)
        
        flowSetUpComplete (1,1){mustBeNonnegative, mustBeNumericOrLogical, mustBeInteger, mustBeLessThanOrEqual(flowSetUpComplete,1)} = false;
        grid; % structures holding gridded interpolants for Cl and Cd
        
    end
    
    % should set a list of properties for models to use/apply
    
    % 3D CL/CD OPTIONS
    % stall delay: default off
    
    % UNSTEADY OPTIONS
    % steady: default on
    % quasi-steady: default off
    % attached unsteady: default off
    % full unsteady (dynamic stall): default off
    
    methods
        
        function [obj] = TidalSim(RunClass, AeroFoilClass)
            %TIDALSIM Constructor method takes "RunClass" containing flow,
            %turbine parameters and operating conditions and "AeroFoilClass"
            %which hold the steady aerodynamic coefficients
            
            obj.Run = RunClass;
            obj.Omega = abs(obj.Run.HubVelocity*obj.Run.TipSpeedRatio/obj.Run.RadialCoord(end));                  % rotational speed of blades (rad/s)
            obj.RotationPeriod = (2*pi)/obj.Omega;
            % set phase shift between blades
            shift = 360/obj.Run.Blades;
            phase = 0:obj.Run.Blades-1;
            obj.Phase = deg2rad(shift.*phase);
            obj.AeroFoil = AeroFoilClass;
            
        end
        
        function [] = setDiscretisation(obj)
            %SETDISCRETISATION Sets the temporal and spatial discretisation
            
            
            % Temporal discretisation
            dt=obj.RotationPeriod/obj.Steps;                   % time step;
            obj.Time = 0:dt:obj.Rotations*obj.RotationPeriod - dt;        % time (s)
            obj.Psi = obj.Omega.*obj.Time;               % temporal azimuthal position (deg)
            
            % Spatial discretisation
            obj.BladeSections = 100;                            % number of blade sections
            obj.RadialCoords = linspace(obj.Run.RadialCoord(1),obj.Run.RadialCoord(end),obj.BladeSections);      % radial coordinates (m)
            obj.BladeTwist = interp1(obj.Run.RadialCoord,obj.Run.BladeTwist,obj.RadialCoords,'PCHIP') + obj.Run.PitchAngle;       % gemoetrical twist + pitch (rad)
            obj.BladeChord = interp1(obj.Run.RadialCoord,obj.Run.BladeChord,obj.RadialCoords,'PCHIP');             % chord length (m)            
            obj.Radius = (obj.RadialCoords(end));                             % radius of blade (m)
            
            
        end
        
        function [] = setOnsetFlow(obj)
            %SETONSETFLOW Calculates the flow due to environmental and operating conditions
            %calls to compute turbiulent fluctuations and wave particle velocity oscillations
            % Components are comnined to form axial and tangential flow
            % seen by the blade
            
            % TODO: PREALLOCATION
            
            % set the flow depending on the operating conditions
            
            % Turbulence
            ut = zeros(size(obj.Time)); wt = ut; vt = ut;
            
            if obj.Run.Turbulence.On
                
                fs=1/(obj.RotationPeriod/obj.Steps); % sampling frequency
                f0 = 1/obj.Time(end);                 % lowest frequency
                f = [f0:f0:fs/2]';                      % frequency array for turbulence model
                df = median(diff(f));                   % if distributed unifoirmly df=f0
                
                % Turbulence is assumed to be perfectly correlated in space
                % (fluctuations at blade1 = bladeN for all N)
                Suvw = vKSpec(obj.Run.HubVelocity,obj.Run.Turbulence.Intensity,obj.Run.Turbulence.LengthScale,...
                    f',obj.Run.Turbulence.IsotropyRatio);% von Karmen PSD for [u,v,w] velocity fluctuations [ m^2/s]
                
                % set or save seed for random number generator             
                if ~isempty(obj.SeedTurbulence)
                    rng(obj.SeedTurbulence)
                    psi = rand(size(Suvw))*2*pi; % random phase
                else
                    psi = rand(size(Suvw))*2*pi; % random phase
                    obj.SeedTurbulence = rng;
                end
                               
                % create velocity fluctuations
                ut = ones(obj.BladeSections,1).*sum(sqrt(2*df*Suvw(:,1)).*cos(2*pi*f.*obj.Time +psi(:,1)));
                vt = ones(obj.BladeSections,1).*sum(sqrt(2*df*Suvw(:,2)).*cos(2*pi*f.*obj.Time +psi(:,2)));
                wt = ones(obj.BladeSections,1).*sum(sqrt(2*df*Suvw(:,3)).*cos(2*pi*f.*obj.Time +psi(:,3)));
            end
            
            % Blade position
            for n = 1:obj.Run.Blades
                z_blade = obj.RadialCoords'.*sin(obj.Psi-obj.Phase(n));
                z = obj.Run.HubDepth + z_blade;
                
                % Shear profile (as sampled by the blade)
                U_shear(n,:,:) = obj.Run.HubVelocity.*((abs(obj.Run.HubDepth) + z_blade)./abs(obj.Run.HubDepth)).^(obj.Run.ShearLaw);
                
                
                if obj.Run.Waves.On
                    wave = Waves; % make waves class
                    wave.WaveCurrent = obj.Run.Waves.WaveCurrent;
                    wave.Type = obj.Run.Waves.Type;
                    wave.Model = obj.Run.Waves.Model;
                    wave.Hs = obj.Run.Waves.Height;
                    wave.Tp = obj.Run.Waves.Period;
                    wave.Direction = obj.Run.Waves.Direction;
                    wave.U0 = obj.Run.HubVelocity;
                    wave.zCord = z;
                    wave.Depth = obj.Run.Depth;
                    
                    if strcmp('Irregular', wave.Type)
                        wave.Periods = obj.Run.Waves.Periods;
                    end
                    
                    % set seed for random number generator   
                    if ~isempty(obj.SeedWaves)
                        wave.Seed = obj.SeedWaves;
                    end
                    
                    % time shift due to wave direction and/or yaw angle
                    tX=obj.RadialCoords'.*sin(obj.Psi-obj.Phase(n))*sin(obj.Run.YawAngle + obj.Run.Waves.Direction)/abs(obj.Run.HubVelocity);      
                    wave.Time = obj.Time + tX;
                    wave.MakeWaves; % run wave model
                    u_wave = wave.UVel; w_wave = wave.WVel; % get the wave partical velocity components
                    obj.Run.Waves.WaveNumber = wave.WaveNumber;
                    obj.SeedWaves = wave.Seed;
                    
                    
%                     [u_wave(n,:,:),w_wave(n,:,:),K] = wavePV(obj.Run.Waves.Height,obj.Run.Waves.Period,obj.Run.Waves.Direction...
%                         ,obj.Run.Depth,z,(obj.Psi-obj.Phase(n)),obj.RadialCoords,obj.Run.TipSpeedRatio,obj.Run.HubVelocity,obj.Run.YawAngle);
%                     obj.Run.Waves.WaveNumber = K;
                else
                    u_wave = zeros(size(U_shear)); w_wave = u_wave;
                end
                
                
                % Total streamwise velocity seen by the blade
                U_axial(n,:,:) = (squeeze(U_shear(n,:,:)) + ut).*cos(obj.Run.YawAngle) ...
                    + squeeze(u_wave(n,:,:)).*cos(obj.Run.YawAngle + obj.Run.Waves.Direction);
                
                % Tangential (DEPTHWISE CONTRIBUTION)
                W_Tan =(squeeze(w_wave(n,:,:)) + wt).*sin(obj.Psi-obj.Phase(n));
                
                % Tangential (TRANSVERSE HORIZONTAL CONTRIBUTION)
                V_Tan =vt.*cos(obj.Psi-obj.Phase(n));
                
                % Yaw and wave angle transverse components
                U_tan =(squeeze(U_shear(n,:,:)) + ut).*sin(obj.Run.YawAngle).*cos(obj.Psi-obj.Phase(n))....
                    + squeeze(u_wave(n,:,:)).*sin(obj.Run.YawAngle + obj.Run.Waves.Direction).*cos(obj.Psi-obj.Phase(n));
                
                % Total tangential velocity
                U_theta(n,:,:) = obj.Omega.*obj.RadialCoords' + U_tan + W_Tan + V_Tan;
            end
            
            obj.UAxial = U_axial;
            obj.UTangential = U_theta;
        end
        
        function [] = setRelativeFlow(obj)
            
            % BEM induction calculation based on steady blade 1 only
            % (uniform current then yaw correction post solution)
            [a,ap] = bladeEM(obj.Run.HubVelocity,obj.Run.TipSpeedRatio,obj.Run.Blades,obj.RadialCoords,....
                obj.BladeChord,obj.BladeTwist,obj.grid);
            
            obj.AxialIndFactor = a;
            obj.TangentialIndFactor = ap;
            
            ap_psi = ap'.*ones(size(obj.Psi));
            a_psi = a'.*ones(size(obj.Psi));
            chi=(0.6.*ap_psi+1)*obj.Run.YawAngle;       % wake skew angle
            Mu = obj.RadialCoords/obj.Radius;
            
            for n = 1:obj.Run.Blades
                % yaw correction to axial induction factor to give azimuthal variation
                ak = a_psi.*(1+((15.*pi)/64).*tan(0.5.*chi).*Mu'.*sin(obj.Psi-obj.Phase(n))); % correction to axial induction factor for skewed wake
                ak(1,:)=0.990;
                ak(end,:)=0.990;
                
                % Relative velocity
                Wrel(n,:,:) = sqrt((squeeze(obj.UAxial(n,:,:)).*(1-ak)).^2 + (squeeze(obj.UTangential(n,:,:)).*(1+ap_psi)).^2);
                
                % Angle of attack
                AoA(n,:,:) = pi/2 - atan2(squeeze(obj.UTangential(n,:,:)).*(1+ap_psi),squeeze(obj.UAxial(n,:,:)).*(1-ak))-obj.BladeTwist';
                
                % Flow angle
                phi(n,:,:) = AoA(n,:,:) + obj.BladeTwist;
                
            end
            
            obj.URelative = Wrel;
            obj.AngleOfAttack = AoA;
            obj.FlowAngle = phi;
            %
        end
        
        function [] = applyRotationalAugmentation(obj)
            
            % Rotational flow augmentation
            F=sepPoint(obj.AeroFoil.AngleOfAttack, obj.AeroFoil.ZeroLiftAngle,obj.AeroFoil.NormCoeff, ...
                obj.AeroFoil.LinearLiftSlope, obj.AeroFoil.LinearClRange); % for rotational solution
            
            [Cl_3d,Cd_3d,~,~] = stallDelay(obj.BladeTwist, obj.RadialCoords, obj.BladeChord, obj.AeroFoil.AngleOfAttack',F' ,...
                obj.AeroFoil.LiftCoeff', obj.AeroFoil.DragCoeff', obj.AeroFoil.LinearLiftSlope, obj.AeroFoil.ZeroLiftAngle);
            
            % Deep stall applied to rotational values (-180 <-> 180)
            for n =1:length(obj.RadialCoords)
                [Values_360r.Alpha, Values_360r.Cl(:,n), Values_360r.Cd(:,n),Values_360r.Cn(:,n)] = vitExtrapolation(obj.AeroFoil.AngleOfAttack,Cl_3d(:,n)',Cd_3d(:,n)');
            end
            % Seperation point for seperation values
            Values_360r.F = sepPoint(Values_360r.Alpha', obj.AeroFoil.ZeroLiftAngle ,Values_360r.Cn, obj.AeroFoil.LinearLiftSlope, obj.AeroFoil.LinearClRange); % Rotational separation point (-180 <-> 180)
            
            % set the rotational coefficients for use
            obj.AeroCoeffs = Values_360r;
            
        end % apply rotationalAugmentation
        
        function [] = setQSBladeCoefficients(obj)                     
            
            % Quasi-steady load analysis: interpolate from gridded tables           
            
            if ~obj.RotationalAugmentation
                % without rotational augmentation
                for n = 1:obj.Run.Blades
                    Cl_QS(n,:,:) = obj.grid.FCl(squeeze(obj.AngleOfAttack(n,:,:)));
                    Cd_QS(n,:,:) = obj.grid.FCd(squeeze(obj.AngleOfAttack(n,:,:)));
                end               
            end
            
            if obj.RotationalAugmentation
                % with rotational augmentation
                outRange=(obj.RadialCoords > 0.8*obj.Radius); % boolean (TRUE for outer section)
                rr = obj.RadialCoords'.*ones(size(squeeze(obj.AngleOfAttack(1,:,:))));
                for n = 1:obj.Run.Blades
                    Cl_QS(n,:,:) = obj.grid.FCl(squeeze(obj.AngleOfAttack(n,:,:)), rr);
                    Cd_QS(n,:,:) = obj.grid.FCd(squeeze(obj.AngleOfAttack(n,:,:)), rr);
                    % Cd near the tip = 2D (No rotatonal augmentation)
                    Cd_QS(n,outRange,:) = obj.grid.Cd2d(squeeze(obj.AngleOfAttack(n,outRange,:)));                                       
                end              
            end
            
            obj.LiftCoeff = Cl_QS; obj.DragCoeff = Cd_QS;
            
                       
        end % Compute quasi-steady blade coefficients
        
        
        function [] = setUSBladeCoefficients(obj)
            
            % set gridded seperation point
            if obj.RotationalAugmentation
                gridF = griddedInterpolant({obj.AeroCoeffs.Alpha, obj.RadialCoords}, obj.AeroCoeffs.F, 'linear');
            else
                gridF = griddedInterpolant(obj.AeroCoeffs.Alpha, obj.AeroCoeffs.F, 'linear');
            end
            
            % pass AoA history to indicial load model
            dt = obj.RotationPeriod/obj.Steps;  % time step;
            for n = 1:obj.Run.Blades
                % Attached unsteady Cl solution
                [Cl_us(n,:,:), Cl_c, Cl_nc, Ds, aE(n,:,:)] = wag(obj.BladeChord, dt, obj.Run.HubVelocity,...
                    squeeze(obj.AngleOfAttack(n,:,:)), obj.AeroFoil.ZeroLiftAngle, obj.AeroFoil.LinearLiftSlope);
                
                if isempty(obj.DSData)
                     error('Dynamic stall empirical data not set. Set the file path on DSData.')
                end
                % Dynamic stall Cl solution
                [~,~,Cl_US(n,:,:), Dvis, Cd_Ind, ff_3d(n,:,:), fff_3d(n,:,:), VortexTracker_3d(n,:,:)] =.....
                    dynStall(obj.BladeTwist, obj.BladeChord, gridF, obj.RadialCoords, squeeze(Cl_us(n,:,:)), ...
                    Cl_c, Cl_nc, Ds, squeeze(aE(n,:,:)), squeeze(obj.AngleOfAttack(n,:,:)), obj.DSData);
                
                % Dynamic stall Cd solution
                [Cd_US(n,:,:)] = dynStallCd(Dvis, Cd_Ind, obj.grid, squeeze(aE(n,:,:)), obj.RadialCoords, obj.AeroFoil.ZeroLiftDrag);
            end
            
            obj.LiftCoeff = Cl_US; obj.DragCoeff = Cd_US;
            
        end
        
        function [] = setRotorLoads(obj, varargin)
            
            [opts, ~] = checkOptions({{'Quasi-steady'},{'Unsteady'}}, varargin);
            
            if ~obj.flowSetUpComplete
                obj.RunSimulation;
            end
            
            if opts(1)
                obj.setQSBladeCoefficients;
                
            elseif opts(2)
                
                obj.setUSBladeCoefficients;
            end
            
            
            for n = 1:obj.Run.Blades
                FF=0.5.*obj.Density.*obj.BladeChord.*obj.URelative(n,:,:).^2; % Dynamic pressure
                [MY(n,:),MX(n,:),T(n,:),P(n,:),FN(n,:,:),FT(n,:,:)]=loads(FF(:,2:end-1,:),...
                    obj.LiftCoeff(n,2:end-1,:),obj.DragCoeff(n,2:end-1,:),obj.FlowAngle(n,2:end-1,:),obj.RadialCoords(2:end-1),obj.Omega);
            end
            
            ext = zeros(obj.Run.Blades,1,size(obj.Psi,2));  % insert zero end values
            obj.ForceNormal = [ext FN ext]; obj.ForceTangential = [ext FT ext];
            obj.RootBM = MY; obj.EdgeBM = MX;
            obj.Power = P; obj.Thrust = T;
            
        end
        
        
        
        function [] = RunSimulation(obj)
            
            obj.setDiscretisation;
            
            % set look-up tables
            if obj.RotationalAugmentation
                % make call here               
                obj.applyRotationalAugmentation;
                obj.grid.FCl = griddedInterpolant({obj.AeroCoeffs.Alpha, obj.RadialCoords}, obj.AeroCoeffs.Cl, 'spline');
                obj.grid.FCd = griddedInterpolant({obj.AeroCoeffs.Alpha, obj.RadialCoords}, obj.AeroCoeffs.Cd, 'spline');
                % need 2d Cl solution for blade sections near the tip
                obj.grid.Cd2d = griddedInterpolant(obj.AeroFoil.CoeffValues_360.Alpha, obj.AeroFoil.CoeffValues_360.Cd, 'spline');
                
            else
                obj.AeroCoeffs = obj.AeroFoil.CoeffValues_360;
                obj.grid.FCl = griddedInterpolant(obj.AeroCoeffs.Alpha, obj.AeroCoeffs.Cl, 'spline');
                obj.grid.FCd = griddedInterpolant(obj.AeroCoeffs.Alpha, obj.AeroCoeffs.Cd, 'spline');
            end
            
            obj.setOnsetFlow;
            obj.setRelativeFlow; % relative velocity, angle of attack, flow angle (BEM called here)
            
            obj.flowSetUpComplete = true; % flag that the flow has been set          
            
            obj.setRotorLoads(obj.LoadMethod);
            
        end % Run simulation
               
    end % methods
    
end