classdef blade < handle
    
    % BLADE In-code representation of a wind turbine blade design
    %   Contains all variables necessary to define the blade geometry
    
    % The point of this (slightly overkill, I admit) object-oriented
    % structure is separation-of-concerns, that is, each file deals with
    % one level of processing only. The blade is not concerned with what
    % the elements do, the elements are not concerned with what the
    % aerofoil lookup tables do, and so on.
    
    % Other things to note:
    % - ALL ANGLES ARE IN RADIANS!
    % (The only exceptions are the lookups inside the aerofoil class and
    % displaying of angles, but those are usually converted on-the-fly)
    
    properties(SetAccess = immutable)
        % Design parameters
        % These are things that we can change about the blade design
        R; % Outside radius of the blade, in m
        tipSpeedRatio;
        rootChord; % Chord length at the tip of the blade, in m
        tipChord; % Chord length at the tip of the blade, in m
        rootAlpha; % Angle of attack at the root of the blade, in RADIANS
        tipAlpha; % Angle of attack at the root of the blade, in RADIANS
        aerofoil; % The aerofoil profile the blade uses
        
        % Derived blade parameters
        rotationSpeed; % Rotation speed of the turbine in rad/s
        discArea; % Area swept out by the turbine in m^2
    end
    
    % Non-final fields
    properties
        solved = false; % True once the blade has been solved
        elements; % Vector of bladeElement objects representing the blade geometry
    end
    
    methods
        
        function this = blade(R, tipSpeedRatio, rootChord, tipChord, rootAlpha, tipAlpha, aerofoil)
            
            % BLADE Creates a new turbine blade with the given parameters
            %   Detailed explanation goes here
            
            this.R = R;
            this.tipSpeedRatio = tipSpeedRatio;
            this.rootChord = rootChord;
            this.tipChord = tipChord;
            this.rootAlpha = rootAlpha;
            this.tipAlpha = tipAlpha;
            this.aerofoil = aerofoil;
            
            chords = linspace(rootChord, tipChord, constants.elementsPerBlade);
            alphas = linspace(rootAlpha, tipAlpha, constants.elementsPerBlade);
            
            r = linspace(constants.hubRadius, this.R, constants.elementsPerBlade);
            
            % First one has to be outside the loop because MATLAB is weird like that
            this.elements = bladeElement(this, r(1), aerofoil, chords(1), alphas(1));
            
            for n = 2:constants.elementsPerBlade
                this.elements(n) = bladeElement(this, r(n), aerofoil, chords(n), alphas(n));
            end
            
            % Derived blade parameters
            this.rotationSpeed = (tipSpeedRatio * constants.windSpeed) / R;
            this.discArea = pi * R^2;
            
        end
        
        function Cp = solve(this)
            
            % SOLVE Solves all elements of this blade design for pitch
            % angle and calculates the resulting output power coefficient
            
            Mtheta = zeros([1, length(this.elements)]);
            
            % Solve individual elements for their moments
            % Vectorised form seems to cause trouble so let's just use a
            % loop for now
            for n = 1:length(this.elements)
                Mtheta(n) = this.elements(n).solve;
            end
            %Mtheta = solve(this.elements);

            % Integrate numerically using the trapezium rule
            T = trapz([this.elements.r], Mtheta); % Square brackets get ALL element radii in one vector!
            
            % Finally, calculate Cp
            Cp = (this.rotationSpeed * T * constants.numBlades) / (0.5 * constants.rho * constants.windSpeed^3 * this.discArea);
            
            % Sanity checks
            if Cp > constants.betzLimit
                error("Rotor Cp %.3f exceeds the Betz limit!", Cp);
            end
            
            if Cp < 0
                error("Rotor Cp %.3f is negative!", Cp);
            end
            
        end
        
        function child = breed(this, that)
            
            % BREED 'Breeds' this blade design with the given blade,
            % combining their parameters and applying random mutations to
            % produce a child blade.
            
            inheritance = rand > 0.5; % Bit of a nasty hack but it's neater than loads of ifs
            childLambda = inheritance .* this.tipSpeedRatio + ~inheritance .* that.tipSpeedRatio;
            if rand > constants.mutationChance
                childLambda = childLambda .* (1 + (rand*2-1) * constants.mutationAmount);
            end
            
            inheritance = rand > 0.5;
            childTipChord = inheritance .* this.tipChord + ~inheritance .* that.tipChord;
            if rand > constants.mutationChance
                childTipChord = childTipChord .* (1 + (rand*2-1) * constants.mutationAmount);
            end
            
            inheritance = rand > 0.5;
            childRootAlpha = inheritance .* this.rootAlpha + ~inheritance .* that.rootAlpha;
            if rand > constants.mutationChance
                childRootAlpha = childRootAlpha .* (1 + (rand*2-1) * constants.mutationAmount);
            end
            
            inheritance = rand > 0.5;
            childTipAlpha = inheritance .* this.tipAlpha + ~inheritance .* that.tipAlpha;
            if rand > constants.mutationChance
                childTipAlpha = childTipAlpha .* (1 + (rand*2-1) * constants.mutationAmount);
            end
            
            % Everything with this. in front of it doesn't change
            child = blade(this.R, childLambda, this.rootChord, childTipChord, childRootAlpha, childTipAlpha, this.aerofoil);
            
        end
        
        function h = visualise(this)
            
            % VISUALISE Plots a 3D view of this blade, showing the
            % individual profiles arranged on the z-axis to give a 3D
            % representation of the shape. Returns a handle to the figure.
            
            h = plot3([0, 0], [0, 0], [constants.hubRadius, this.R], 'b');
            
            hold on;
            grid on;
            axis equal;
            view(250, 30);
            xlim([-0.02, constants.maxChord + 0.02]);
            ylim([-0.02, constants.maxChord + 0.02]);
            zlim([constants.hubRadius, this.R]);
            
            for n = 1:length(this.elements)
                this.elements(n).plot3D;
            end
                
        end
        
        function [] = export(this, folderpath)
            
            % EXPORT Exports this blade as a set of CSV files in the given
            % folder. Each element will have its own separate file of xyz
            % coordinates describing the blade profile at that position.
            
            if ~isfolder(folderpath)
                error("Folder %s does not exist, please create it first!");
            end
            
            for n = 1:length(this.elements)
                this.elements(n).export(strcat(folderpath, "/Element_", num2str(n), ".csv"));
            end
            
        end
        
    end
    
end

