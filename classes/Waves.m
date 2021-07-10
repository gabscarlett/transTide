classdef Waves < handle
    % WAVES: Computes wave properties; spectrum, linear dispersion 
    % and wave partical velocities for regular or irregular waves
    % to first or second order. 

    properties (Access = public)
        Type(1,1) string {mustBeMember(Type, ["Regular","Irregular"])} = 'Regular'; % wave type
        Model(1,1) string {mustBeMember(Model, ["Linear","Second"])} = 'Linear'; % wave model (linear wave theory or Stokes second order)
        Hs(1,1){mustBeNumeric, mustBeFinite}; % significant wave height
        Tp(1,1){mustBeNumeric, mustBeFinite}; % peak wave period
        Direction (1,1){mustBeNumeric, mustBeFinite}; % wave direction   
        Periods (1,:){mustBeNumeric, mustBeFinite} = []; % relative wave period array for irregular waves
        Spectra(1,1) string {mustBeMember(Spectra, ["Bretschnider","Jonswap"])} = 'Bretschnider'; % wave type;
        Depth (1,1){mustBeNumeric, mustBeFinite}= 45; % water depth
        zCord(:,:){mustBeNumeric, mustBeFinite} = []; % depth position array
        Time(:,:){mustBeNumeric, mustBeFinite} = []; % simulation time array
        U0(1,1){mustBeNumeric, mustBeFinite}= 2; % streamwise current velocity;
        W0(1,1){mustBeNumeric, mustBeFinite}= 0; % depthwise current velocity;;
        
    end
    
    properties(SetAccess = private)
        
        UVel;
        WVel;
        Spectrum;
        Amplitudes;        
        WaveNumber;
             
    end
    
    properties (Access = private)        
        modOpt = 0; % used to switch wave equation between 1st and 2nd order        
    end
    
    
    methods
        function [obj] = MakeWaves(obj)
            
            if strcmp('Second', obj.Model)
                obj.modOpt = 1; % switches on 2nd order component
            end
            
            
            if strcmp('Regular', obj.Type)
                obj.computeVelocity(obj.Tp, obj.Hs);
            end
            
            if strcmp('Irregular', obj.Type)
                obj.makeSpectrum;
                phase = (2*pi*rand(1, length(obj.Periods)));
                obj.computeVelocity(obj.Periods, 2*obj.Amplitudes, 'Phase', phase);
            end
            
            
        end
        
        function [obj] = makeSpectrum(obj)
            
            w_r = 2*pi./obj.Periods;
            dw = w_r(2:end) - w_r(1:end-1);
            dw(end+1) = dw(end);   
            
            wp = 2*pi./obj.Tp;
            coef = 5/16*wp^4*obj.Hs^2;          
            obj.Spectrum = coef./(w_r.^5).*exp(-5/4*(wp./w_r).^4);
            
            obj.Amplitudes = sqrt(2*obj.Spectrum.*abs(dw));%*length(obj.S);
            
        end
        
        function [obj] = computeVelocity(obj, T, H, varargin)
            [opts, args] = checkOptions({{'Phase',1}}, varargin);
            
            phase = 0;
            g =9.81; x = 0;
            w_r = 2*pi./T;
            
            if opts(1)
                phase = args{1};
            end
            
            % call wavenumber function to iterate for k
            tol=1E-10;     % iteration error
            k = wavenumber(g, w_r,obj.Direction, obj.Depth, obj.U0, obj.W0, tol);
            T_a=2*pi./sqrt(k.*g.*tanh(k.*obj.Depth));
            w_a = 2*pi./T_a; % doppler shifted frequency
            
            t_wave = obj.Time;
            
            for n = 1:length(H)
                U(:,:,n) =(H(n)*g*k(n)/(2*w_a(n))).*(cosh(k(n)*(obj.Depth + obj.zCord))/cosh(k(n)*obj.Depth)).*cos(k(n)*x-w_a(n).*t_wave + phase(n)) +....
                    obj.modOpt*(3*H(n)^2*w_a(n)*k(n)/16).*(cosh(2*k(n).*(obj.Depth + obj.zCord))./(sinh(k(n)*obj.Depth))^4).*cos(2*(k(n)*x-w_a(n).*t_wave + phase(n)));
                
                W(:,:,n) =(H(n)*g*k(n)/(2*w_a(n))).*(sinh(k(n)*(obj.Depth + obj.zCord))/cosh(k(n)*obj.Depth)).*sin(k(n)*x-w_a(n).*t_wave + phase(n)) +....
                    obj.modOpt*(3*H(n)^2*w_a(n)*k(n)/16).*(sinh(2*k(n).*(obj.Depth + obj.zCord))./(sinh(k(n)*obj.Depth))^4).*sin(2*(k(n)*x-w_a(n).*t_wave + phase(n)));
            end
            obj.UVel = squeeze(sum(U,3));
            obj.WVel = squeeze(sum(W,3));
            obj.WaveNumber = k;
            
        end
        
    end
    
end