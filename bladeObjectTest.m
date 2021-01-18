% Test script for the object-oriented blade solver structure

clear;
close all;

% Design parameters
R = 0.5; % Turbine radius (m)
tipSpeedRatio = 2.5; % Sensible magnitude 1-10
rootAlpha = deg2rad(5); % Desired root angle of attack (radians)
tipAlpha = deg2rad(10); % Desired tip angle of attack (radians)
rootChord = 0.2; % Root chord length (m)
tipChord = 0.1; % Tip chord length (m)

% Initialise the aerofoil object
S822 = aerofoil("S822", "aerofoil_data");

tic;

% Initialise the blade
b = blade(R, tipSpeedRatio, rootChord, tipChord, rootAlpha, tipAlpha, S822);

% Solve and display Cp
Cp = b.solve;
fprintf("Cp = %.3f\n", Cp);

fprintf("Solver run time = %.3fs\n", toc); % Run time readout

% Display the blade in 3D
b.visualise;