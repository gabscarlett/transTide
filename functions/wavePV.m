function [u_wave,w_wave,k]=wavePV(Hs,Tw,wd,h,z,PSI,r,TSR,U0,gamma)

% FUNCTION THAT DETERMINES WAVE INDUCED VELOCITIES AT EACH BLADE LOCATION


%% INPUTs

    
    % Hs=5.0;      % significant wave height (m)
    % Tw=10.0;     % zero crossing period (s)
    % h=46;        % water depth, free surface to seabed (m) 
    % PSI =(omega.*t - phase); 
    
    %wd=0;                  % direction of wave (rad)
    x=0;                    % [0:lambda/100:lambda] horizontal coordinate (m);
    W0=0;                   % depthwise current (m/s)
    g=9.81;                 % acceleration due to gravity (m/s^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigma_r=2*pi/Tw; % angular frequency of wave (rad/s)

    % call wavenumber function to iterate for k
    tol=1E-10;     % iteration error
    k = wavenumber(g,sigma_r,wd,h,U0,W0,tol);
    
    Ta=2*pi/sqrt(k*g*tanh(k*h));
    sigma = 2*pi/Ta; % doppler shifted frequency

        
        R=r(end);
        % yaw correction to wave model
        tI=(PSI)*R/(TSR*abs(U0));                   % time at blade          
        tX=r'.*sin(PSI)*sin(gamma+wd)/abs(U0);      % time shift due to wave direction and/or yaw angle
        t_wave=tI+tX;                               % total wave time
        

        % wave particle velocities at time, and depth z - Stokes second-order wave theory
        % with yaw correction
        u_wave =(Hs*g*k/(2*sigma)).*(cosh(k*(h+z))/cosh(k*h)).*cos(k*x-sigma.*t_wave)+(3*Hs^2*sigma*k/16).*(cosh(2*k.*(h+z))./(sinh(k*h))^4).*cos(2*(k*x-sigma.*t_wave));
        w_wave =(Hs*g*k/(2*sigma)).*(sinh(k*(h+z))/cosh(k*h)).*sin(k*x-sigma.*t_wave)+(3*Hs^2*sigma*k/16).*(sinh(2*k.*(h+z))./(sinh(k*h))^4).*sin(2*(k*x-sigma.*t_wave));
        
        
end