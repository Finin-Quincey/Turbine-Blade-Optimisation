% Various DoEs for the turbine blade (factorial, latin hypercube)

clear;
close all;

% Initialise the aerofoil object
S822 = aerofoil("S822", "aerofoil_data");

% Design space to explore
minR = 0.5;
maxR = 0.5;
minLambda = 2.9153;
maxLambda = 2.9153;
minAlpha = deg2rad(8.5);
maxAlpha = deg2rad(9.5);
minChord = 0.095;
maxChord = 0.105;
%radii = linspace(0.3, 0.5, 5); % Turbine radius (m)
%tipSpeedRatios = linspace(0.5, 5, 5); % Sensible magnitude 1-10
%alphas = deg2rad(linspace(1, 10, 5)); % Desired angles of attack (radians)
%chords = linspace(0.05, 0.25, 5); % Chord lengths (m)

% Initialise parameter indices for the full-factorial DoE runs
%runs = fullfact(5, 5, 5, 5); % Unlike lhsdesign, this supports different numbers of levels
    
% Normalisation
% Full factorial outputs *levels* rather than fractions, so we need this
% extra step to get it in the same format as lhsdesign (I get why, but for
% this it's irrelevant)
% Divide each column by the max value in that column
%for k = 1:size(runs, 2)
%    runs(:, 1) = runs(:, 1) / max(runs(:, 1));
%end

runs = lhsdesign(1000, 4, 'iterations', 20); % A thousand runs takes around 40 seconds

N = size(runs, 1); % Total number of runs
X = zeros(N, 4); % Init input matrix for use in analysis later
Cp = zeros(N, 1); % Init output vector for use in analysis later

% Profiling
runtimes = zeros(N, 1); % Init runtime graph data
startTime = tic; % Start overall timer

% Profiler plot setup
runtimeGraph = bar(runtimes, 'FaceColor', 'flat', 'EdgeColor', 'none');
profiler = gcf;
grid on;
xlim([0, N]);
xlabel("Run #");
ylabel("Run duration (s)");

%figure;
%bladeVis = axes;

for n = 1:N % For each run of the DoE
    
    % Calculate parameter values from normalised runs
    R = minR + (maxR - minR) * runs(n, 1);
    lambda = minLambda + (maxLambda - minLambda) * runs(n, 2);
    alpha = minAlpha + (maxAlpha - minAlpha) * runs(n, 3);
    chord = minChord + (maxChord - minChord) * runs(n, 4);
    
    X(n, :) = [R, lambda, alpha, chord];
    
    runStart = tic; % Start individual run timer (excludes the above setup stuff)
    
    % I'm just using constant chord and alpha for now
    b = blade(R, lambda, 0.2, chord, alpha, alpha, S822);
    
    try
        Cp(n) = b.solve; % Try to solve the blade for this set of values
        %figure(bladeVis);
        %b.visualise;
    catch
        Cp(n) = NaN; % Use NaN for all errors (no convergence, mathematical invalidity, out-of-bounds Re/alpha, etc.)
    end
    
    % Update profiler
    
    runtimes(n) = toc(runStart);
    
    pause(0.001); % Let the graph sort itself out
    
    if ~isvalid(runtimeGraph)
        error("Terminated by user"); % Close the profiler window to abort
    end
    
    runtimeGraph.YData = runtimes;
    if isnan(Cp(n))
        runtimeGraph.CData(n, :) = [0.8, 0, 0]; % Colour failed runs red...
    else
        runtimeGraph.CData(n, :) = [0.1, 0.8, 0]; % ... and successful runs green
        %cla(bladeVis); % Super slow, super cool
    end
    
    title(sprintf("Elapsed time: %.2fs", toc(startTime)));
    
end

close(profiler);

% Display failed run percentage (helps figure out how effectively we're
% using our resources!)
fprintf("Proportion of failed runs: %.1f%%\n", sum(isnan(Cp)) / N * 100);

interactionAnalysis(X, Cp);