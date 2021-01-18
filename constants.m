classdef constants
    
    % CONSTANTS Stores a bunch of constants for the turbine blade
    % optimisation, so they can be accessed from anywhere.
    
    properties(Constant = true)
        
        % Design constants
        windSpeed = 3; % Wind speed in m/s
        numBlades = 3; % Number of times I'm going to make them print my blade
        hubRadius = 0.12; % Inside radius of the blade in m
        elementsPerBlade = 30; % Number of elements along the blade length
        maxChord = 0.2; % Maximum allowed chord length in m
        
        % Physical constants
        rho = 1.225; % Density of air at 20degC, 1bar (kg/m^3)
        mu = 1.8205e-5; % Dynamic viscosity of air at 20degC, 1bar
        betzLimit = 0.593; % Maximum possible value of Cp
        
        % Solver constants
        convergenceThreshold = 0.001; % Convergence threshold as a fraction of the value
        iterationLimit = 50; % Number of iterations before solver aborts
        
        % GA constants
        mutationChance = 0.5; % Chance of a given blade parameter being mutated
        mutationAmount = 0.02; % Fractional change in values when mutated
        
    end
    
end

