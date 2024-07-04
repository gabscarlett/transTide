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

%% make an AerofoilProps class by passing the file name
foil = AerofoilProps(fileNameFoil);

%% make a RunConditions class by passing the turbine file and operating file
run = RunConditions('turbine file',fileNameTurb, 'operating file', fileNameOps);

% to turn off turbulence/waves;
%run.Turbulence.On = false;
%run.Waves.On = false;

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
sim.Rotations = 5; % reduce the number of simulated rotations from 100 to 50
%sim.Steps = 36; % reduce the steps in each rotation from 72 (5 degrees) to 36 (10 degrees)

%% run a simulation using default settings

sim.RunSimulation; % run the simulation

Power_0 = sum(sim.Power); % sum the power contribution of each plade
meanPower_0 = mean(Power_0); % compute the mean power

Thrust_0 = sum(sim.Thrust); % sum the thrust contribution of each plade
meanThrust_0 = mean(Thrust_0); % compute the mean thrust

meanThrust_0/1000; % print the mean thrust in kN

%% apply pitch control and rerun

% This example inputs a sinusoidal dynamic pitch. This is in addition to
% the static pitch.

% get the time array
t = sim.Time; % get the time array

% Set the ampitude of the signal
pitchAmp = 0.5; % amplitude in degrees

% frequency of the signal
f = 1/sim.RotationPeriod; % this is the rotational frequency,

% loop over the blades (we're only shifting the phase)
for n = 1: length(sim.Phase)
    pitchDyn(n,:) = deg2rad(pitchAmp) * sin(2*pi*f.*t + sim.Phase(n)); % dynamic pitch time series in radians
end

% now set the pitch control on the simulation.

sim.PitchControl.On = true;
sim.PitchControl.Pitch = pitchDyn;

%% rerun the simulation and compare power an thrust

sim.RunSimulation; % run the simulation

Power_p = sum(sim.Power); % sum the power contribution of each plade
meanPower_p = mean(Power_p); % compute the mean power

Thrust_p = sum(sim.Thrust); % sum the thrust contribution of each plade
meanThrust_p = mean(Thrust_p); % compute the mean thrust

meanThrust_p/1000; % print the mean thrust in kN

%% plot the power and thrust time series

figure(1);
% plot power
subplot(1,2,1)
plot(t, Power_0/1000, 'k', t, meanPower_0/1000*(ones(size(t))), 'k:', 'LineWidth',2)
hold on
plot(t, Power_p/1000, 'r', t, meanPower_p/1000*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Power [kW]')
% plot thrust
subplot(1,2,2)
plot(t, Thrust_0/1000, 'k', t, meanThrust_0/1000*(ones(size(t))), 'k:', 'LineWidth',2)
hold on
plot(t, Thrust_p/1000, 'r', t, meanThrust_p/1000*(ones(size(t))), 'r:', 'LineWidth',2)
xlabel('Time [s]')
ylabel('Thrust [kN]')
%print([myPath 'plots\' 'power_and_thrust_pitchControl'],'-dpng','-r500');
