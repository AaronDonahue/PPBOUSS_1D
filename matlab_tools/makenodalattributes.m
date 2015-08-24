%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Script to construct nodal attributes file				%
%										%
% This Matlab code has been designed to quickly construct a nodal attributes    %
% file. This controls things such as the bottom friction, sponge generation and %
% sponge absorption.								%
%										%
% The following is a list of the component that must be changed for each new	%
% attributes file.								%
%										%
% nodefile : Name of the nodal attributes file, usually *****.13		%
% noderun : Name for the run, this is written on the first line of the file.    %
% p : This is the order of the run, usually set to p = 1, linear.		%
% gridfile : The grid file to be used to construct the attributes file.  	%
%										%
% Note, in order for this to work one must go to the "Nodal Attributes" section %
% and change key terms, these are the;						%
% attr(ind).on : Either 0 or 1, if the attribute will be used then this should  %
%                be set to 1, otherwise 0 means don't use.			%
% attr(ind).default : This is the default nodal value.				%
%										%
% Furthermore, one needs to adjust terms in each of the subsections associated  %
% with each individual attribute, such as sponge generation or absorbtion.      %
% This gives the user full control over each attribute.  At the end of the run  %
% plots for each attribute are plotted, so that the user can determine if the   %
% attribute is set in the way that is desired. 					%
%										%
% There are a few pre-loaded hotstart setups which can be uncommented.		%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nodal attributes file control variables
% nodefile = 'fort.13';
% noderun = 'test';
% p = 1;
% gridfile = 'fort.14';

% --Ting and Kirby
%  nodefile = 'fort.13';
%  noderun = 'TingKirby';
%  p = 1;
%  gridfile = 'TingKirby.14';

% --Solitary Wave (breaking and nonbreaking)
%  nodefile = 'SolitaryWave.13'
%  noderun = 'Solitary Wave'
%  p = 1;
%  gridfile = 'Solitary.14';

% --Breakwater (Hsiao)
 nodefile = 'Hsiao.13';
 noderun = 'Hsiao';
 p = 1;
 gridfile = 'Hsiao_case1.14';
%  %  gridfile = 'Hsiao_case2.14';
%  %  gridfile = 'Hsiao_case3.14';

% --Carrier and Greenspan
%  nodefile = 'carrier.13';
%  noderun = 'Carrier_Greenspan';
%  p = 1;
%  gridfile = 'carrier.14';

% --Dambreak
%  nodefile = 'dambreak.13';
%  noderun = 'dambreak';
%  p = 1;
%  gridfile = 'dambreak.14';

% --Riemann
%  nodefile = 'riemann.13';
%  noderun = 'riemann';
%  p = 1;
%  gridfile = 'riemann.14';

% --Parabolic
%  nodefile = 'parabolic.13';
%  noderun = 'parabolic';
%  p = 1;
%  gridfile = 'parabolic.14';

%% Read grid file
fid = fopen(gridfile);
tmp = textscan(fid,'%f %f',1,'headerlines',1);
ne = tmp{1};
nn = tmp{2};
tmp = textscan(fid,'%f %f %f',nn);

xt = tmp{2};
dt = tmp{3};

x = zeros(ne,p+1);
for l = 1:ne
    x(l,:) = linspace(xt(l),xt(l+1),p+1);
end

%% Nodal Attributes
ind = 0;
% Bottom Friction
ind = ind+1;
attr(ind).name      = 'mannings_n_at_sea_floor';
attr(ind).unit      = 'm^2/s';
attr(ind).on        = 1;
attr(ind).default   = 0.004;
attr(ind).numnondef = 0;
attr(ind).nondefvalue = [];
attr(ind).nondefelem  = [];

% Sponge Generation
ind = ind+1;
attr(ind).name    = 'sponge_generation_layer';
attr(ind).unit    = 'unitless';
attr(ind).on      = 0;
attr(ind).default = 0;
attr(ind).numnondef = 0;
attr(ind).nondefvalue = [];
attr(ind).nondefelem  = [];

% Sponge Absorbing
ind = ind+1;
attr(ind).name    = 'sponge_absorbing_layer';
attr(ind).unit    = 'unitless';
attr(ind).on      = 1;
attr(ind).default = 0;
attr(ind).numnondef = 0;
attr(ind).nondefvalue = [];
attr(ind).nondefelem  = [];

%
totalattr = ind;
%
nwp = 0;
for i = 1:totalattr
  nwp = nwp + attr(i).on;
end

%% Sponge Generation
if attr(2).on == 1
    samp = 30;
    slen = 10;
    sord = 3;
    x0 = xt(1);
    sfun = @(x)samp/slen*(sord+1)*(1-(x-x0)/slen).^sord;
    for l = 1:ne
        if xt(l) < x0+slen
            attr(2).numnondef = attr(2).numnondef + 1;
            attr(2).nondefvalue = [attr(2).nondefvalue;sfun(x(l,:))];
            attr(2).nondefelem = [attr(2).nondefelem;l];
        end
    end
end

%% Sponge Absorbing
if attr(3).on == 1
    samp = 30;
    slen = 5;
    sord = 3;
    x0 = xt(1);
    sfun = @(x)samp/slen*(sord+1)*(1-(x-x0)/slen).^sord;
%     x0 = xt(end);
%     sfun = @(x)samp/slen*(sord+1)*((x-x0+slen)/slen).^sord;
    for l = 1:ne
        if xt(l) < x0+slen %xt(l) > x0-slen
            attr(3).numnondef = attr(3).numnondef + 1;
            attr(3).nondefvalue = [attr(3).nondefvalue;sfun(x(l,:))];
            attr(3).nondefelem = [attr(3).nondefelem;l];
        end
    end
end



%% Write attributes file
fid  = fopen(nodefile,'w');
fprintf(fid,'%s\n',noderun);
fprintf(fid,'%d %d\n',ne,nn);
fprintf(fid,'%d\n',nwp);

for i = 1:totalattr
    if attr(i).on == 1
        fprintf(fid,'%s\n',attr(i).name);
        fprintf(fid,'%s\n',attr(i).unit);
        fprintf(fid,'%d\n',p+1);
        for j = 1:p+1
          fprintf(fid,'%f ',attr(i).default);
        end
        fprintf(fid,'\n');
    end
end

for i = 1:totalattr
    if attr(i).on == 1
        fprintf(fid,'%s\n',attr(i).name);
        fprintf(fid,'%d\n',attr(i).numnondef);
        for j = 1:attr(i).numnondef
          fprintf(fid,'%d ',attr(i).nondefelem(j));
          for k = 1:p+1
              fprintf(fid,'%f ',attr(i).nondefvalue(j,k));
          end
          fprintf(fid,'\n');
        end        
    end
end


fclose(fid);

%% Visualize Solution
for i = 1:totalattr
    if attr(i).on == 1
        figure(i)
        plot(x(attr(i).nondefelem,:),attr(i).nondefvalue)
        title(attr(i).name)
        xlim([xt(1) xt(end)])
    end
end