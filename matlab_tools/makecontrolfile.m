%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Matlab script to construct a control file                       %
%                                                                         %
% This Matlab code has been designed to quickly construct a control file  %
% file for use by the DG_WASUPP pressure-Poisson based model.             %
%                                                                         % 
% To use this script simply go through each subsection and change the     %
% appropriate value to reflect the value you desire.                      %
%                                                                         %
% Updated: Aug 21, 2015                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Control Parameters
% Run Name (This is a descriptor of the simulation)
runname = 'test';

% Grid File (The name of the grid file to be used with the simulation)
% *REQUIRED
gridfile = 'fort.14';

% Hotstart File (These two variables control whether or not a hotstart is
% going to be used.
ihot = 0;  % 0 = No Hot Start, 1 = Hot Start
hotstart_file = 'fort.67'; % If a hotstart is going to be used this is the 
                           % filename, required if ihot = 1

% Output Directory (change this if you want to output to a different
% directory)
out_direc = './';

% Boundary type (This is the boundary condition flag)
% Currently "reflective" and "radiant" are supported.
boundarytype = 'reflective';

% Run parameters *REQUIRED
p = 1;             % Order of method (integer)
nrk = 3;           % Runge-Kutta method (3 is RKSSP 3,3 for example)
negp = 8;          % Number of Gauss integration points to use, (less means faster, could lose accuracy)
dgbasis = 'nodal'; % Type of DG run, default to nodal, usually stick with this.

% Time control variable *REQUIRED
maxtime = 10;    % Maximum time for run
cfl = 0.5;       % CFL condition, max should be 1 for stability
timesnap = -999; % Time snap for global output, negative means no output

% Wetting/Drying parameters
iwet = 0;     % Wetting Drying off/on, 0 = off, 1 = on
h0 = 0.00001; % Minimum water depth for wetting and drying, required if iwet = 1

% Slope limiter
islp = 0;          % Slope limiter off/on, 0 = off, 1 = on
islpconstant = 50; % "M" value in slope limiter, determines a minimum threshold
                   % for when the slope limiter will be used, this is
                   % neccessary to avoid spurious flattening of true peaks
                   % and troughs.

% Breaking Model
ibreak = 0;  % Breaking on/off, 0 = off, 1 = on
breakingmodel = 'duran'; % Breaking model to be used, currently
                         % tonelli = Depth based, 0.8 of depth
                         % duran = DG based, measures edge discontinuity
                         % duran_adj = Like duran, but measure internal
                         %             slope of element to determine
                         %             trouble elements.

% Nodal Attributes (i.e. friction, sponge generation/absorption)
nwp = 2; % Number of nodal attributes to be used, must match nodal attributes file
nodalattr_file = 'fort.13'; % Name of nodal attributes file. Required if nwp>0

% Station output
station_file = 'none'; % Name of file containing stations.  This file follows
                               % a simple format.  It is as follows:
                               % L1        : Number of stations (integer - N)
                               % L2 - LN+1 : Station location in terms of x
station_timestep = -999;       % Time snap for station output;
                               % If value is less than dt then will snap
                               %   every timestep.

                               
% Pressure Poisson control
inonhydro = 2; % Level of pressure-Poisson solver, 0 = SWE, 2 = mu^2, 4 = mu^4


%% Write control file
fid = fopen('fort.wasupp','w');

fprintf(fid,'!---------------------------------------------------!\n');
fprintf(fid,'! Control file for %s simulation\n',runname);
fprintf(fid,'!---------------------------------------------------!\n');
fprintf(fid,'! Grid and Hotstart Files\n');
fprintf(fid,'grid_file = %s\n',gridfile);
fprintf(fid,'ihot = %d\n',ihot);
fprintf(fid,'hotstart_file = %s\n',hotstart_file);
fprintf(fid,'out_direc = %s\n',out_direc);
fprintf(fid,'!---------------------------------------------------!\n');
fprintf(fid,'! DG Control Variables\n');
fprintf(fid,'p = %d\n',p);
fprintf(fid,'dgbasis = %s\n',dgbasis);
fprintf(fid,'negp = %d\n',negp);
fprintf(fid,'nrk = %d\n',nrk);
fprintf(fid,'boundarytype = %s\n',boundarytype);
fprintf(fid,'!---------------------------------------------------!\n');
fprintf(fid,'! Time control variables\n');
fprintf(fid,'maxtime = %.16f\n',maxtime);
fprintf(fid,'timesnap = %.16f\n',timesnap);
fprintf(fid,'cfl = %.16f\n',cfl);
fprintf(fid,'!---------------------------------------------------!\n');
fprintf(fid,'! Station file\n');
fprintf(fid,'station_file = %s\n',station_file);
fprintf(fid,'station_timestep = %.16f\n',station_timestep);
fprintf(fid,'!---------------------------------------------------!\n');
fprintf(fid,'! Nodal Attributes\n');
fprintf(fid,'nwp = %d\n',nwp);
fprintf(fid,'nodalattr_file = %s\n',nodalattr_file);
fprintf(fid,'!---------------------------------------------------!\n');
fprintf(fid,'! Physics \n');
fprintf(fid,'iwet = %d\n',iwet);
fprintf(fid,'h0 = %.16f\n',h0);
fprintf(fid,'islp = %d\n',islp);
fprintf(fid,'islpconstant = %.16f\n',islpconstant);
fprintf(fid,'ibreak = %d\n',ibreak);
fprintf(fid,'breakingmodel = %s\n',breakingmodel);
fprintf(fid,'!---------------------------------------------------!\n');
fprintf(fid,'! pressure-Poisson\n');
fprintf(fid,'inonhydro = %d\n',inonhydro);
fprintf(fid,'!---------------------------------------------------!\n');

fclose(fid);



