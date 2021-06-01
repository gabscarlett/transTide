clear; close all; clc

% Example script to demonstrate how to run a tidal simulation from file.
% Input data for the turbine, aerofoil and operating conditions are loaded from
% from a csv file.

% Default discretisation settings are changed.

% An example is given to compute and plot the power, thrust and root and edge
% wise bending moments.

%% set the paths and file names
myPath = [transTidePath '\data\SupGenTankScale\']; % this is the data path

dataNameTurb = 'turbine_SupGen'; 
fileNameTurb = [myPath dataNameTurb]; % details of blade profile

dataNameOps = 'operationalConditions_SupGen';
fileNameOps = [myPath dataNameOps]; % details of operating conditions (flow and turbine)

dataNameFoil = 'static_aerofoil_NACA_63_816';
fileNameFoil = [myPath dataNameFoil]; % measured aerofoil coefficients with angle of attack

dsData = 'NACA_63_816_DS_parameters.mat'; % empirical data for dynamic stall model (valid for S814 only)
dsFile = [myPath dsData]; % file path for DS data

%% make an AerofoilProps class by passing the file name
foil = AerofoilProps(fileNameFoil);

%% make a RunConditions class by passing the turbine file and operating file
run = RunConditions('turbine file',fileNameTurb, 'operating file', fileNameOps);

%% adjust flow settings

%run.Turbulence.On = 0; % switch off turbulence
%run.Waves.On = 0; % switch off waves

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

sim.BladeSections = 20; % for the tank scale device 20 sections are plenty
sim.Rotations = 50; % reduce the number of simulated rotations from 100 to 50

%% run a simulation using default settings

sim.RunSimulation; % run the simulation

blade = 1; % blade number to inspect

RootBM = sim.RootBM(blade,:); % root bending moment time series for blade 
meanRootBM = mean(RootBM); % compute the mean root bending moment

EdgeBM = sim.EdgeBM(blade,:); % root bending moment time series for blade 
meanEdgeBM = mean(EdgeBM); % compute the mean root bending moment

Power = sum(sim.Power); % sum the power contribution of each plade
meanPower = mean(Power); % compute the mean power

Thrust = sum(sim.Thrust); % sum the thrust contribution of each plade
meanThrust = mean(Thrust); % compute the mean thrust


%% plot the bending moment time series

t = sim.Time; % get the time array

figure;
% plot power
subplot(2,2,1)
plot(t, RootBM, 'k', t, meanRootBM*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Root bending moment [Nm]')
% plot thrust
subplot(1,2,2)
plot(t, EdgeBM, 'k', t, meanEdgeBM*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Edge wise bending moment [Nm]')

%% plot power and thrust time series

t = sim.Time; % get the time array

figure;
% plot power
subplot(1,2,1)
plot(t, Power, 'k', t, meanPower*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Power [W]')
% plot thrust
subplot(1,2,2)
plot(t, Thrust, 'k', t, meanThrust*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Thrust [N]')