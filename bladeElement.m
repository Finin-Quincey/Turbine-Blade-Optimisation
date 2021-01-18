classdef bladeElement < handle
    
    % BLADEELEMENT Represents a single radial element of a turbine blade.
        
    % Final fields
    properties(SetAccess = immutable)
        parentBlade; % A reference to the parent blade object this element belongs to
        r; % The radial position of this element, in m
        aerofoil; % Aerofoil object
        chord; % Chord length in m
        alpha; % Angle of attack in RADIANS
    end
    
    % Non-final fields
    properties
        solved = false; % True once the profile has been solved
        points; % TRANSFORMED aerofoil points
        pitch; % Pitch angle in RADIANS
    end
    
    methods
        
        function this = bladeElement(parentBlade, r, aerofoil, chord, alpha)
            
            % BLADEELEMENT Creates a new blade element with the given parameters.
            %
            %   bladeElementparentBlade, r, aerofoil, chord, alpha) creates
            %   an element with the given parent blade object, radial
            %   position, aerofoil geometry, pitch angle and chord length.
            %
            %   parentBlade is a reference to the parent blade object this
            %   element belongs to
            %
            %   r is the radial position of this element, in m
            %
            %   aerofoil is a reference to the aerofoil object that
            %   represents the blade profile to be used
            %
            %   alpha is the desired angle of attack, in RADIANS
            %
            %   chord is the chord length of the profile, in m
            
            this.parentBlade = parentBlade;
            this.r = r;
            this.aerofoil = aerofoil;
            this.chord = chord;
            this.alpha = alpha;
        end
        
        function Mtheta = solve(this)
            
            % SOLVE Solves this element for the rotational force component
            % Ftheta. This also sets the pitch angle.
            
            % Starting values
            axInd = 1/3; % Axial induction a
            angInd = 0; % Angular induction a'
            inflowAngle = 0;
            
            iterations = 0;
            converged = false;
            
            blade = this.parentBlade; % Assign to local var for conciseness
            
            while ~converged
    
                converged = true; % It has converged unless we find it hasn't

                % Calculate some dimensionless thingies
                localSpeedRatio = blade.tipSpeedRatio * this.r / blade.R;
                localBladeSolidity = (constants.numBlades * this.chord) / (2 * pi * this.r);

                % Store existing value so we can find the difference later
                prevInflowAngle = inflowAngle;
                % 1. Calculate inflow angle from a and a'
                inflowAngle = atan((1 - axInd) / (localSpeedRatio * (1 + angInd)));

                % Can't divide by zero, and the convergence test is
                % meaningless first time round anyway
                if prevInflowAngle ~= 0
                    % Convergence test
                    converged = converged && abs(prevInflowAngle - inflowAngle) / inflowAngle < constants.convergenceThreshold;
                end

                % Calculate relative wind and hence Reynold's number
                Vrel = (constants.windSpeed * (1 - axInd)) / sin(inflowAngle);
                Re = (constants.rho * Vrel * this.chord) / constants.mu;

                % Look up CL and CD
                % Our nice shiny aerofoil object will take care of this for us
                CL = this.aerofoil.lookupCL(Re, this.alpha);
                CD = this.aerofoil.lookupCD(Re, this.alpha);

                % Be careful which way round the cos and sin are!
                Cx = CL * cos(inflowAngle) + CD * sin(inflowAngle);
                Ctheta = CL * sin(inflowAngle) - CD * cos(inflowAngle);

                % Store existing values so we can find the difference later
                prevAxInd = axInd;
                prevAngInd = angInd;

                % Calculate the next a and a'
                axInd = 1 / ((4 * sin(inflowAngle) ^ 2) / (localBladeSolidity * Cx) + 1);
                angInd = 1 / ((4 * sin(inflowAngle) * cos(inflowAngle)) / (localBladeSolidity * Ctheta) - 1);

                % More sanity checks
                if axInd > 0.5
                    error("Axial induction > 0.5, this is not good!");
                end

                % More convergence tests
                converged = converged && abs(prevAxInd - axInd) / axInd < constants.convergenceThreshold;
                converged = converged && abs(prevAngInd - angInd) / angInd < constants.convergenceThreshold;

                iterations = iterations + 1; % Next iteration
                
                if iterations > constants.iterationLimit % Failsafe
                    error("Blade element at radius %.2fm failed to converge!", this.r);
                end
                    
            end
            
            this.solved = true; % Yay we did it
            
            % Assign the newly-calculated values to the object's fields
            
            % Calculate pitch from inflow angle to give desired alpha
            this.pitch = inflowAngle - this.alpha;
            
            Vrel = (constants.windSpeed * (1 - axInd)) / sin(inflowAngle); % Do we need this?
            Mtheta = 0.5 * constants.rho * Vrel^2 * this.chord * this.r * Ctheta; % Easier to multiply by r in here
            
            % Finally, assign the transformed points
            
            % Scale the profile
            pts = [this.aerofoil.points(:, 1) * this.chord, this.aerofoil.points(:, 2) * this.chord];
            % Now rotate it
            pts = [pts(:, 1) * sin(this.pitch) + pts(:, 2) * cos(this.pitch), pts(:, 2) * sin(this.pitch) - pts(:, 1) * cos(this.pitch)];
            % Finally flip it because we want it to go clockwise
            this.points = [pts(:, 1), -pts(:, 2)];
            
        end
        
        function [] = print(this)
            
            % PRINT Prints the parameters of this blade element to the
            % console in readable form.
            
            fprintf("Blade element at radius %.2fm:\n", this.r);
            
            fprintf("==============================\n");
            
            fprintf("Aerofoil profile: %s\n", this.aerofoil.name);
            fprintf("Angle of attack: %.1f degrees\n", rad2deg(this.alpha));
            fprintf("Chord length: %.2fm\n", this.chord);
            
            if this.solved
                
                fprintf("Solved\n");
            
                fprintf("==============================\n");

                fprintf("Pitch angle: %.1f degrees\n", rad2deg(this.pitch));
                
            else
                fprintf("Not yet solved\n");
            end
            
        end
        
        function h = visualise(this)
            
            % VISUALISE Plots a 2D view of this profile, scaled and rotated
            % appropriately. Returns a handle to the figure.
            
            if ~this.solved
                error('Profile must be solved for pitch angle first!');
            end
            
            h = plot(this.points(:, 1), this.points(:, 2), 'k'); % Profile outline
            hold on;
            scatter(this.points(:, 1), this.points(:, 2), 'b.'); % Profile points
            plot([0, sin(this.pitch) * this.chord], [0, cos(this.pitch) * this.chord], 'c'); % Chord line
            xlabel('x (m)');
            ylabel('y (m)');
            grid on;
            axis equal;
        end
        
        function [] = plot3D(this)
            
            % PLOT3D Plots the geometry of this blade element onto the
            % current set of 3D axes. Used as part of the 3D blade
            % visualisation.
            
            if ~this.solved
                error('Profile must be solved for pitch angle first!');
            end
            
            plot3(this.points(:, 1), this.points(:, 2), ones(size(this.points, 1), 1) * this.r, 'k'); % Profile outline
            plot3([0, sin(this.pitch) * this.chord], [0, cos(this.pitch) * this.chord], [this.r, this.r], 'c'); % Chord line
            
        end
        
        function [] = export(this, filepath)
            
            % EXPORT Exports this blade element to a CSV file containing
            % xyz coordinates describing the outline of the blade profile
            % at this element's position.
            
            % Don't save the first and last points, they can cause
            % geometry issues
            % Multiply by 100 to get it in cm, seems to be what the Fusion
            % plugin uses
            writematrix([this.points(2 : end-1, :), ones(size(this.points, 1) - 2, 1) * this.r] * 100, filepath);
            
        end
        
    end
end

