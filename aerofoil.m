classdef aerofoil % No need for this to be a handle class
    
    % AEROFOIL Represents an aerofoil profile. Deals with loading the data,
    % lift/drag coefficient lookup and storing geometry.
    
    % Yeah okay I'm only ever going to have one of these but it saves
    % having unwieldy lookup tables cluttering up the solver classes
    
    properties(SetAccess = immutable) % No changey
        name; % Name of this profile
        ReValues; % Reynold's number vector (column headers)
        alphaValues; % Angle of attack vector (row headers)
        CLLookup; % Lift coefficient lookup table
        CDLookup; % Drag coefficient lookup table
        points; % 2D coordinates describing the outline of this profile
    end
    
    methods
        
        function this = aerofoil(name, folder)
            
            % AEROFOIL Creates a new aerofoil with the given name and reads
            % the relevant data from the given directory. folder should be
            % a (relative) path to a directory containing .mat files for CL
            % lookup, Re values, alpha values and coordinate points.
            
            this.name = name;
            
            stem = strcat(folder, '/', name);
            
            % Loading .mat files is painful and I'm lazy so this will do
            this.alphaValues = readFirst(strcat(stem, '_Alpha_Values.mat'));
            this.ReValues = readFirst(strcat(stem, '_Re_Values.mat'));
            this.CLLookup = readFirst(strcat(stem, '_CL_Lookup.mat'));
            this.CDLookup = readFirst(strcat(stem, '_CD_Lookup.mat'));
            
            % You just had to be different, didn't you?
            pts = readFirst(strcat(stem, '_Points.mat'));
            this.points = [pts.x, pts.y];
            
            % Prepare the grid for interp2
            [this.ReValues, this.alphaValues] = meshgrid(this.ReValues, this.alphaValues);
            
        end
        
        function CL = lookupCL(this, Re, alpha)
            
            % LOOKUPCL Returns the lift coefficient for the given Reynold's
            % number and angle of attack (in RADIANS) for this aerofoil
            
            % Sanity checks (and boy, do I need my sanity checked...)
            if Re < min(this.ReValues, [], 'all') || Re > max(this.ReValues, [], 'all')
                error("Reynold's number %i is out of bounds!", Re);
            end
            
            if alpha < min(this.alphaValues, [], 'all') || alpha > max(this.alphaValues, [], 'all')
                error("Angle of attack %.1f degrees is out of bounds!", rad2deg(alpha));
            end
            
            % Convert alpha to degrees because the lookup tables are in
            % degrees, not radians!
            CL = interp2(this.ReValues, this.alphaValues, this.CLLookup, Re, rad2deg(alpha));
            
        end
        
        function CD = lookupCD(this, Re, alpha)
            
            % LOOKUPCD Returns the drag coefficient for the given Reynold's
            % number and angle of attack (in RADIANS) for this aerofoil
            
            % Sanity checks
            if Re < min(this.ReValues, [], 'all') || Re > max(this.ReValues, [], 'all')
                error("Reynold's number %i is out of bounds!", Re);
            end
            
            if alpha < min(this.alphaValues, [], 'all') || alpha > max(this.alphaValues, [], 'all')
                error("Angle of attack %.1f degrees is out of bounds!", rad2deg(alpha));
            end
            
            %CD = 0; % Previously we assumed no drag
            CD = interp2(this.ReValues, this.alphaValues, this.CDLookup, Re, rad2deg(alpha));
            
        end
        
    end
    
end

