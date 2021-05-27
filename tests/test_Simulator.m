% test preprocessor

clear; close all; clc

myPath = [transTidePath '\data\'];
dataNameTurb = 'turbine_TGL';
dataNameOps = 'operationalConditions';
fileNameTurb = [myPath dataNameTurb];
fileNameOps = [myPath dataNameOps];

dataNameFoil = 'static_aerofoil_S814';
fileNameFoil = [myPath dataNameFoil];

foil = AerofoilProps(fileNameFoil);

run = RunConditions('turbine file',fileNameTurb, 'operating file', fileNameOps);

sim = TidalSim(run, foil); % pass the run settings class and aerofol class to simulator class

sim.LoadMethod = "Unsteady";

%sim.RotationalAugmentation = true;

dsData = 'S814_DS_parameters.mat';

dsFile = [myPath dsData];

sim.DSData = dsFile;

sim.RunSimulation;