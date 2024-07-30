clear; close all; clc

% Example script to demonstrate how to run a tidal simulation from file.
% Input data for the turbine, aerofoil and operating conditions are loaded from
% from a csv file.

% An example is given to make a pitch time series and to set it on the
% class to control the pitch angle

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

% to turn off turbulence/waves;
run.Turbulence.On = false;
run.Waves.On = false;

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

%% how to change discretisation settings

sim.BladeSections = 50; % reduce blade dicretisation from 100 to 50 sections
sim.Rotations = 1; % reduce the number of simulated rotations from 100 to 50
%sim.Steps = 36; % reduce the steps in each rotation from 72 (5 degrees) to 36 (10 degrees)

%% run a simulation using default settings

sim.RunSimulation; % run the simulation

%% loop over static pitch angles

% The static pitch is set using Run.PitchAngle, or in the excel input file.

pitchStatic = 0:0.1:5;

for k = 1: length(pitchStatic)

    run.PitchAngle = ones(1,length(sim.Phase))*deg2rad(pitchStatic(k));

    sim.RunSimulation; % run the simulation

    pow = sum(sim.Power); % sum the power contribution of each blade
    powMean(k) = mean(pow);
    powStd(k) = std(pow);

end

%% plot the mean and standard deviation with pitch amplitude

subplot(1,2,1)
plot(pitchStatic, powMean, 'k', 'LineWidth',2)
xlabel('pitch amplitude [deg]')
ylabel('average power [W]')
grid on
subplot(1,2,2)
plot(pitchStatic, powStd, 'k', 'LineWidth',2)
xlabel('pitch amplitude [deg]')
ylabel('standard deviation of power [W]')
grid on
