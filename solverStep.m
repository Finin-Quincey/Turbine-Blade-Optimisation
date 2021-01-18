function T = solverStep(rMin, R, N, B, c, axInd, angInd, tipSpeedRatio)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Ftheta = zeros(N, 1);

n = 1;

% Loop through each element of the span
for r = linspace(rMin, R, N)
    
    localSpeedRatio = tipSpeedRatio * r/R;
    
    localBladeSolidity = (B * c) / (2 * pi * r);
    
    inflowAngle = atan((1 - axInd) / (localSpeedRatio * (1 + angInd)));
    
%     pitchAngle = inflowAngle - angleOfAttack;
    
    Vrel = (windSpeed * (1 - axInd)) / sin(inflowAngle);
    
    rei = (rho * Vrel * c) / mu;
    
    % Look up CL and CD
    % BE CAREFUL, THIS FILE IS IN DEGREES, NOT RADIANS!!!
    CL = interp2(Re, alpha, S822_CL_Lookup, rei, rad2deg(angleOfAttack));
    CD = 0; % For now assume CD = 0, will change once we have the data
    
    Cx = CL * sin(inflowAngle) + CD * cos(inflowAngle);
    Ctheta = CL * cos(inflowAngle) - CD * sin(inflowAngle);
    
    axInd = 1 / ((4 * sin(inflowAngle) ^ 2) / (localBladeSolidity * Cx) + 1);
    angInd = 1 / ((4 * sin(inflowAngle) * cos(inflowAngle)) / (localBladeSolidity * Ctheta) - 1);
    
    Ftheta(n) = 0.5 * rho * Vrel^2 * c * r * Ctheta;
    
    n = n + 1;
    
end

T = trapz(Ftheta);

end

