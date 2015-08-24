%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Matlab script to construct a grid file                                %
%                                                                               %
% This Matlab code has been designed to quickly construct a grid                %
% file for use by the DG_WASUPP pressure-Poisson based model.                   %
%                                                                               % 
% The following is a list of components that must be changed   			%
% for each new grid file: 							%
% 										%
% runname : The is a string that will be written on the first line 		%
%	    of the grid file to identify the file. 				%
% gridfile : This is the actual name of the file, usually ****.14. 		%
% dx : The element size, currently setup for structured meshes, so 		%
%      only one dx value is needed. 						%
% xl,xr : The left and right endpoints for the mesh, respectively. 		%	 
%   										%
% Farther down in the code one also needs to define the bathymetry.   		%
% this is done by changing the 'dep_fun' variable as a function of 'x'.  	%
% For example dep_fun = @(x)0*x + 1 will make a flat bathymetry with depth 1.   %
%                                                                               %
% There are a few pre-loaded mesh setups which can be uncommented to get those. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Controls for gridfile
%  runname = 'test';
%  gridfile = 'fort.14';
%  dx = 1;
%  xl = -10;
%  xr = 10;
%  dep_fun = @(x)0*x+1;


% --Ting and Kirby
% runname = 'Ting_Kirby';
% gridfile = 'TingKirby.14';
% dx = 0.025;
% xl = -25;
% xr = 25;
%  hb = 0.4;
%  hb_0 = 0;
%  hb_slope = 1/35;
% dep_fun = @(x)0*x+hb - .02/.2*(x-hb_0).*(x>=hb_0).*...
%           (x<hb_0+0.2)-(0.02+hb_slope*(x-hb_0-0.2)).*(x>=hb_0+0.2);

% --Solitary Wave (Synolakis)
% runname = 'SolitaryWave_break';
% gridfile = 'Solitary.14';
%  dx = 0.25;
%  xl = -25;
%  xr = 25;
%  dep_fun = @(x)1 - 1/19.85*(x+19.85).*(x>=-19.85);

% --Breakwater (Hsiao)
runname = 'Hsiao_case1';
gridfile = 'Hsiao_case1.14';
 dx = 0.02;
 xl = -10;
 xr = 25;
hb = 0.18;
eadj = 0.02;
% eadj = 0.04;
% eadj = 0.076;
dep_fun = @(x) (hb - 1/20*(x-10).*(x>=10).*(x<=13.6) ...
       - (1/4*(x-13.6)+hb).*(x>13.6).*(x<13.9) ...
       - (0.076+hb).*(x>=13.9).*(x<=13.948) ...
       - (-1/1.8*(x-13.948)+0.076+hb).*(x>=13.948).*(x<=14.045) ...
       - 1/20*(x-10).*(x>14.045))+eadj;

% --Carrier and Greenspan problem
% runname = 'Carrier_Greenspan';
% gridfile = 'carrier.14';
%  dx = 0.1;
%  xl = -45;
%  xr = 5;
% dep_fun = @(x)1 - 1/19.85*(x+19.85).*(x>=-19.85);

% --Dambreak problem
% runname = 'dambreak';
% gridfile = 'dambreak.14';
%  dx = 10;
%  xl = -300;
%  xr = 300;
% dep_fun = @(x)0*x;

% --Riemamnn Problem
% runname = 'riemann';
% gridfile = 'riemann.14';
%  dx = 1.25
%  xl = -200;
%  xr = 400;
% dep_fun = @(x)0*x;

% --Parabolic Basin
% runname = 'parabolic';
% gridfile = 'parabolic.14';
%  dx = 12.5;
%  xl = -4000;
%  xr = 4000;
% dep_fun = @(x)-1.6*10^(-7)*x.^2;

%% Determine x vector
x = xl:dx:xr;
nn = length(x);

%% Determine depth
dep = dep_fun(x);

% for i = 1:10
%     dep = smooth(dep);
% end

%% Write file
fid = fopen(gridfile,'w');
fprintf(fid,'%s\n',runname);
fprintf(fid,'%d %d\n',nn-1,nn);
for i = 1:nn
    fprintf(fid,'%4d %16.8f %16.8f\n',[i,x(i),dep(i)]);
end
fclose(fid);

%% Plot grid file to check
plot(x,-dep)
