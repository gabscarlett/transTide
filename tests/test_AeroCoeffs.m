% test AerofoilCoefficients script

clear; close all; clc

a = AerofoilCoefficients;

a.Units = 'deg';

a.AngleOfAttack= 1:10;

a.NormCoeff = (1:10) * 0.9 + 0.665;

a.convertCoeffs('Cl')

CL = a.LiftCoeff;
AOA = a.AngleOfAttack;

plot(AOA,CL)