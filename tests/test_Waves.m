% script to test waves

clear all; close all; clc,

%% make wave object for regular wave and compute velocities 

% by default the Class "Waves" is set to run regular waves
wave = Waves;
wave.Hs = 2; % for regular this is relative height
wave.Tp =10; % for regular this is relaive period
wave.U0 = 2; % current velocity
wave.Depth = 45; % water depth

wave.Time = 0:150;
z = -10; % the depthwise coordinate can be a single value in time or SEE LINE BELOW ...
%z = linspace(-10,-20,length(wave.Time)); % vector of vertical positions in time which should be the same length as time.
wave.zCord = z; % position below the free surface to analyse

% compute the wave particle velocities and plot
wave.MakeWaves;

figure;
plot(wave.Time, wave.UVel,'b')
hold on
plot(wave.Time, wave.WVel, 'r')
xlabel('Time [s]')
ylabel('Velocity [m/s]')

%% make wave object for irregular wave and compute velocities

wave = Waves;
% change property "WaveType" to Irregular
wave.Type = 'Irregular';
wave.Hs = 4; % significant wave height
wave.Tp = 10; % peak wave period
wave.U0 = 2; % current velocity
wave.Periods = 1:0.1:15; % the wave period to examine
wave.Depth = 45; % water depth

% w=2*pi./wave.Periods
% [S, A] = bretschneider(wave.Hs, wave.Tp, w, [])

wave.Time = 0:0.1:1500;
%z = -10; % the depthwise coordinate can be a single value in time or SEE LINE BELOW ...
z = linspace(-10,-20,length(wave.Time)); % vector of vertical positions in time which should be the same length as time.
wave.zCord = z; % position below the free surface to analyse

% compute the wave particle velocities and plot
wave.MakeWaves;

figure;
plot(wave.Time, wave.UVel,'b')
hold on
plot(wave.Time, wave.WVel, 'r')
xlabel('Time [s]')
ylabel('Velocity [m/s]')