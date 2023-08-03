clear; close all; clc

% Example script to demonstrate how to run a tidal simulation from file.
% Input data for the turbine, aerofoil and operating conditions are loaded from
% from a csv file.

% Altering default discretisation settings is demonstrated.

% An example is given to compute the root and edge wise bending moments,
% and to plot them for a single blade.

%% set the paths and file names
myPath = [transTidePath '\data\TGLFullScale\']; % this is the data path

dataNameTurb = 'turbine_TGL'; 
fileNameTurb = [myPath dataNameTurb]; % details of blade profile

dataNameOps = 'operationalConditions_TGL';
fileNameOps = [myPath dataNameOps]; % details of operating conditions (flow and turbine)

dataNameFoil = 'static_aerofoil_S814';
fileNameFoil = [myPath dataNameFoil]; % measured aerofoil coefficients with angle of attack

dsData = 'S814_DS_parameters.mat'; % empirical data for dynamic stall model (valid for S814 only)
dsFile = [myPath dsData]; % file path for DS data

%% make an AerofoilProps class by passing the file name
foil = AerofoilProps(fileNameFoil);

%% make a RunConditions class by passing the turbine file and operating file
run = RunConditions('turbine file',fileNameTurb, 'operating file', fileNameOps);

%% make a TidalSim class by passing the AerofoilProps and Runconditions class

sim = TidalSim(run, foil); % pass the run settings class and aerofol class to simulator class

% TidalSim is the simulation class. It pulls everything together and runs the simulation. 
% Model options and discretisation are controlled within.
% Functions and methods to compute the flow field and loads
% on the tubine are called.

%%%%%%%% Default simulation options %%%%%%%%%%%%%%%%%%

% quasi-steady
% non-rotational

%%%%%%%% Default discretisation %%%%%%%%%%%%%%%%%%

% Number of blade sections: 100
% Number of rotations made: 100
% Steps per revolution: 72 (step size of 5 degrees)

%% change discretisation settings

sim.BladeSections = 50; % reduce blade dicretisation from 100 to 50 sections
sim.Rotations = 50; % reduce the number of simulated rotations from 100 to 50
sim.Steps = 36; % reduce the steps in each rotation from 72 (5 degrees) to 36 (10 degrees)

%% run a simulation using default settings

sim.RunSimulation; % run the simulation

blade = 1; % blade number to inspect
RootBM = sim.RootBM(blade,:); % root bending moment time series for blade 
meanRootBM = mean(RootBM); % compute the mean root bending moment
EdgeBM = sim.EdgeBM(blade,:); % root bending moment time series for blade 
meanEdgeBM = mean(EdgeBM); % compute the mean root bending moment

%% plot the bending moment time series

t = sim.Time; % get the time array

figure;
% plot power
subplot(1,2,1)
plot(t, RootBM/1000, 'k', t, meanRootBM/1000*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Root bending moment [kNm]')
% plot thrust
subplot(1,2,2)
plot(t, EdgeBM/1000, 'k', t, meanEdgeBM/1000*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Edge wise bending moment [kNm]')

%% inspect the bending moments along the blade span

blade = 1; % blade number to inspect

r = sim.RadialCoords;
FN = sim.ForceNormal(blade,:,:);
FT = sim.ForceTangential(blade,:,:);

for n = 2:length(r)-1 % integrate along the blade in increments (dr) from the root to the tip

    MY(n,:)=simps(r(1:n),FN(:,1:n,:).*(r(1:n)-r(1)));                % root bending (N m)
    MX(n,:)=simps(r(1:n),FT(:,1:n,:).*(r(1:n)-r(1)));                % edgewise bending (N m)

end

% plot the summarry statistics