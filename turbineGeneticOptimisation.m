%% Genetic algorithm for optimising the turbine blade

% A spot of housekeeping...
clear;
close all;

%% Setup

% GA parameters
numGenerations = 400; % How many iterations of the genetic algorithm to perform
populationSize = 32; % Number of blades in each generation (should be a multiple of 4)

% Solver parameters
S822 = aerofoil("S822", "aerofoil_data");
R = 0.5;

% Graphics
grey = [0.5, 0.5, 0.5];
blue = [0, 0, 1];
cyan = [0, 1, 1];

saveFrames = true; % True to save images of each generation to a folder

% Starting values
startTipChord = 0.1;
startRootAlpha = deg2rad(4);
startTipAlpha = deg2rad(4);
startLambda = 5;

% Don't ask. MATLAB is weird.
dummyBlade = blade(0, 0, 0, 0, 0, 0, S822);

% Init with dummy blade
population = dummyBlade;

% Initialise population
for n = 1:populationSize
    % Randomise parameters by the global mutation amount
    tipChord = startTipChord * (1 + (rand*2-1) * constants.mutationAmount);
    rootAlpha = startRootAlpha * (1 + (rand*2-1) * constants.mutationAmount);
    tipAlpha = startTipAlpha * (1 + (rand*2-1) * constants.mutationAmount);
    lambda = startLambda * (1 + (rand*2-1) * constants.mutationAmount);
    % Create the resulting blade object
    population(n) = blade(R, lambda, constants.maxChord, tipChord, rootAlpha, tipAlpha, S822);
end

bestCp = -inf; % Will be used to track the best Cp value found

%% Visualisation Setup

% Initialise visualisation
fig = figure('Position', [200, 60, 1200, 700]);
tiledlayout(4, 2);

bladeVisAxes = nexttile([4, 1]);

var1Axes = nexttile;
var1Plot = scatter([], [], '.');
xlabel('Tip Chord Length (m)');
ylabel('C_p');
grid on;
xlim([0, constants.maxChord]);
ylim([0, constants.betzLimit]);

var2Axes = nexttile;
var2Plot = scatter([], [], '.');
xlabel('Root Angle of Attack (deg)');
ylabel('C_p');
grid on;
xlim([0, 12]);
ylim([0, constants.betzLimit]);

var3Axes = nexttile;
var3Plot = scatter([], [], '.');
xlabel('Tip Angle of Attack (deg)');
ylabel('C_p');
grid on;
xlim([0, 12]);
ylim([0, constants.betzLimit]);

var4Axes = nexttile;
var4Plot = scatter([], [], '.');
xlabel('Tip Speed Ratio');
ylabel('C_p');
grid on;
xlim([1, 7]);
ylim([0, constants.betzLimit]);

%% Main Loop

for gen = 1:numGenerations
    
    %% Solver
    
    Cps = zeros(length(population), 1); % Init Cp array for this generation
   
    % Solve each blade in the population for Cp
    for n = 1:length(population)
        try
            Cps(n) = population(n).solve;
        catch
            Cps(n) = -inf; % Failed runs get the lowest possible Cp
        end
    end
    
    %% Update Parameter Visualisations
    
    [~, bestCpIndex] = max(Cps); % Need this before we sort Cp for visualisation
    
    drawnow; % Let the graph sort itself out
    
    if ~isvalid(fig)
        error("Terminated by user"); % Close the window to abort
    end
    
    % Update parameter plots
    var1Plot.XData = [var1Plot.XData, [population.tipChord]];
    var1Plot.YData = [var1Plot.YData, Cps'];
    var1Plot.CData = repmat(grey, populationSize * gen, 1);
    var1Plot.CData(populationSize * (gen-1) + 1 : end, :) = repmat(blue, populationSize, 1);
    var1Plot.CData(populationSize * (gen-1) + bestCpIndex, :) = cyan;
    var1Axes.Title.String = sprintf("Generation %i / %i", gen, numGenerations);
    
    var2Plot.XData = [var2Plot.XData, rad2deg([population.rootAlpha])];
    var2Plot.YData = [var2Plot.YData, Cps'];
    var2Plot.CData = repmat(grey, populationSize * gen, 1);
    var2Plot.CData(populationSize * (gen-1) + 1 : end, :) = repmat(blue, populationSize, 1);
    var2Plot.CData(populationSize * (gen-1) + bestCpIndex, :) = cyan;
    
    var3Plot.XData = [var3Plot.XData, rad2deg([population.tipAlpha])];
    var3Plot.YData = [var3Plot.YData, Cps'];
    var3Plot.CData = repmat(grey, populationSize * gen, 1);
    var3Plot.CData(populationSize * (gen-1) + 1 : end, :) = repmat(blue, populationSize, 1);
    var3Plot.CData(populationSize * (gen-1) + bestCpIndex, :) = cyan;
    
    var4Plot.XData = [var4Plot.XData, [population.tipSpeedRatio]];
    var4Plot.YData = [var4Plot.YData, Cps'];
    var4Plot.CData = repmat(grey, populationSize * gen, 1);
    var4Plot.CData(populationSize * (gen-1) + 1 : end, :) = repmat(blue, populationSize, 1);
    var4Plot.CData(populationSize * (gen-1) + bestCpIndex, :) = cyan;
    
    %% Evolutionary Selection
    
    [Cps, idx] = sort(Cps); % Get the order of indices from smallest Cp to largest
    
    % Eliminate the worst designs by taking only the top half
    population = population(idx(length(idx)/2 + 1:end));
    
    % Track the best blade so far and its corresponding Cp value
    if Cps(end) > bestCp
        
        bestBlade = population(end);
        bestCp = Cps(end);
    
        %% Update Blade Visualisation
        % No point doing this unless it changed!
    
        if ~isvalid(fig)
            error("Terminated by user"); % Close the window to abort
        end

        axes(bladeVisAxes); % Make blade vis the current axes
        cla(bladeVisAxes); % Clear the old blade
        title(sprintf("Current Best C_p = %.4f", bestCp));
        bestBlade.visualise;

        drawnow; % Let the graph sort itself out

    end
    
    % Save animation frames
    if saveFrames
        saveas(fig, strcat("animation_frames/generation_", num2str(gen), ".png"));
    end
    
    %% Next Generation
    
    % Shuffle the order so they breed randomly
    population = population(randperm(length(population)));
    
    % Init with dummy blade
    nextGeneration = dummyBlade;
    
    for n = 1:2:length(population) % Go through the remaining blades in pairs
        % Breed each pair of blades 4 times to keep population size constant
        for k = 1:4
            nextGeneration(2*n - 1 + (k-1)) = population(n).breed(population(n+1));
        end
    end
    
    population = nextGeneration; % Assign the new generation of blades
    
end