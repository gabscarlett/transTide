clear; close all; clc

% Example script to demonstrate how to run a tidal simulation from file.
% Input data for the turbine, aerofoil and operating conditions are loaded from
% from a csv file.

% An example is given to compute the power with and without the effect of
% rotational augmentation (3D stall delay) accounted for.

%% set the paths and file names
myPath = [transTidePath '\data\']; % this is the data path

dataNameTurb = 'turbine_TGL'; 
fileNameTurb = [myPath dataNameTurb]; % details of blade profile

dataNameOps = 'operationalConditions';
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

%% run a simulation using default settings

sim.RunSimulation; % run the simulation

Power = sum(sim.Power); % sum the power contribution of each plade
meanPower = mean(Power); % compute the mean power

Thrust = sum(sim.Thrust); % sum the thrust contribution of each plade
meanThrust = mean(Thrust); % compute the mean thrust


%% re-run the simulation with rotational augmentation

sim.RotationalAugmentation = true; % set rotational augmentation property to true
sim.RunSimulation; % re-run the simulation

Power_R = sum(sim.Power); % sum the power contribution of each blade
meanPower_R = mean(Power_R); % compute the mean power

Thrust_R = sum(sim.Thrust); % sum the thrust contribution of each blade
meanThrust_R = mean(Thrust_R); % compute the mean thrust


%% plot the power and thrust time series for each method

t = sim.Time; % get the time array

figure;
% plot power
subplot(1,2,1)
plot(t, Power/1000, 'b', t, meanPower/1000*(ones(size(t))), 'b:', ...
    t, Power_R/1000, 'r', t, meanPower_R/1000*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Power [kW]')
legend('Without rotation', 'mean','With rotation', 'mean', 'location', 'best')
% plot thrust
subplot(1,2,2)
plot(t, Thrust/1000, 'b', t, meanThrust/1000*(ones(size(t))), 'b:', ...
    t, Thrust_R/1000, 'r', t, meanThrust_R/1000*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Thrust [kN]')