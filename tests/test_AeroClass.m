% test AerofoilCoefficients script

clear; close all; clc

myPath = [myMatPath '\transTide\data\'];
dataName = 'static_aerofoil_S814';
fileName = [myPath dataName];
%foilTab = readtable([fileName '.csv']);

a = AerofoilProps(fileName);
a
