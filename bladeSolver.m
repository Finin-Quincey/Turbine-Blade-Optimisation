% Initial solution script for the main wind turbine blade assignment, used
% to get the solution right first before doing all the fancy stuff (it's
% intentionally a bit messy)

clear;
close all;

debug = false;

% Constants
N = 30; % Number of elements in blade span
R = 0.5; % Turbine radius (m)
B = 3; % Number of blades
rMin = 0.12; % Hub diameter (m)
angleOfAttack = deg2rad(linspace(5, 5, N)); % Desired angles of attack (radians)
c = linspace(0.2, 0.2, N); % Chord length (m)
tipSpeedRatio = 2; % Sensible magnitude 1-10
windSpeed = 3; % Wind speed (m/s)
rho = 1.225;%1.2047; % Density of air at 20degC, 1bar (kg/m^3)
mu = 1.8205e-5; % Dynamic viscosity of air at 20degC, 1bar
% https://theengineeringmindset.com/properties-of-air-at-atmospheric-pressure/
convergenceThreshold = 0.0001; % Convergence threshold as a fraction of the value

% Derived constants
rotationSpeed = (tipSpeedRatio * windSpeed) / R;
discArea = pi * R^2;

% Load aerofoil data
load('aerofoil_data/S822_Alpha_Values.mat');
load('aerofoil_data/S822_Re_Values.mat');
load('aerofoil_data/S822_CL_Lookup.mat');

[Re, alpha] = meshgrid(S822_Re_Values, S822_Alpha_Values);

% Starting values
axInd = 1/3; % Axial induction a
angInd = 0; % Angular induction a'
inflowAngle = 0;

if debug
    Fthetas = zeros(1, N);
end

converged = false([1, N]);
iterations = 0;

% Loop through each element of the span
r = linspace(rMin, R, N);

while sum(converged) ~= length(converged) && iterations < 1000
    
    converged = true([1, N]); % It has converged unless we find it hasn't

    localSpeedRatio = tipSpeedRatio .* r/R;

    localBladeSolidity = (B * c) ./ (2 * pi .* r);

    prevInflowAngle = inflowAngle;
    inflowAngle = atan((1 - axInd) ./ (localSpeedRatio .* (1 + angInd)));

    % Dirty hack for the first element
    if prevInflowAngle ~= 0
        converged = converged & abs(prevInflowAngle - inflowAngle) ./ inflowAngle < convergenceThreshold;
    end

    pitchAngle = inflowAngle - angleOfAttack;

    Vrel = (windSpeed .* (1 - axInd)) ./ sin(inflowAngle);

    rei = (rho .* Vrel .* c) ./ mu;

%     if min(rei) < min(Re) || max(rei) > max(Re)
%         error("Reynolds number out of bounds!");
%     end
    
    % Look up CL and CD
    % BE CAREFUL, THIS FILE IS IN DEGREES, NOT RADIANS!!!
    CL = interp2(Re, alpha, S822_CL_Lookup, rei, rad2deg(angleOfAttack));
    CD = 0; % For now assume CD = 0, will change once we have the data

    % Be careful which way round the cos and sin are!
    Cx = CL .* cos(inflowAngle) + CD .* sin(inflowAngle);
    Ctheta = CL .* sin(inflowAngle) - CD .* cos(inflowAngle);

    prevAxInd = axInd;
    prevAngInd = angInd;

    axInd = 1 ./ ((4 .* sin(inflowAngle) .^ 2) ./ (localBladeSolidity .* Cx) + 1);
    angInd = 1 ./ ((4 .* sin(inflowAngle) .* cos(inflowAngle)) ./ (localBladeSolidity .* Ctheta) - 1);

    if sum(axInd > 0.5) > 0
        error("Axial induction > 0.5, maths is invalid!");
    end
    
    converged = converged & abs(prevAxInd - axInd) ./ axInd < convergenceThreshold;
    converged = converged & abs(prevAngInd - angInd) ./ angInd < convergenceThreshold;

    Vrel = (windSpeed .* (1 - axInd)) ./ sin(inflowAngle);

    Ftheta = 0.5 .* rho .* Vrel.^2 .* c .* r .* Ctheta;
    
    if debug
        Fthetas = [Fthetas; Ftheta];
    end
    
    iterations = iterations + 1;
    
end
    
T = trapz(r, Ftheta);

if sum(converged) ~= length(converged)
    Cp = NaN;
    error("Failed to converge!");
else
    Cp = (rotationSpeed * T * B) / (0.5 * rho * windSpeed^3 * discArea);
end

if debug
    pause(1);
    animatedPlot(@plot, iterations+1, 1:iterations+1, Fthetas, r);
end