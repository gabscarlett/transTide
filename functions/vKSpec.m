function [Suvw] = vKSpec(U,I,L,f,R)

if ~isempty(f(f<0)), error('"frequency" must be positive'); end
% ---------------------------------------------
% INPUTS

% U:        mean (streamwise velocity) [m/s]
% I:        Turbulence intensity in [u,v,w]
% L:     	Length scale
% f:        Frequency array [Hz]
% R:        Ratio Ivw/Iu; 1 = isotropic turbulence
% ---------------------------------------------
% OUTPUTS 

%Suvw:      PSD for [u,v,w] velocity fluctuations [ m^2/s]
% ---------------------------------------------
U=abs(U);
tau =I*U.*ones(size(f)); % Standard deviation
% Von Karman coefficent
coef=[-5/6, -11/6];
% dimension of output
Suvw = zeros(size(f,2),3);
%calculation of S/std^2
nx = L.*U.^(-1).*f;
ny = R*L.*U.^(-1).*f;
% streamwise
Suvw(:,1) =  U.^(-1).*4.*L.*tau.^2.*(1+70.7.*nx.^2).^(coef(1));
% vertical and transverse
Suvw(:,2)=  U.^(-1).*4.*R*L.*(R*tau).^2.*(1+70.7.*4.*ny.^2).^(coef(2)).*(1+188.4.*4*ny.^2);
Suvw(:,3)=Suvw(:,2);

end