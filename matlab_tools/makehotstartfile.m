%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Script to construct hotstart file					%
%										%
% This Matlab code has been designed to quickly construct a hotstart file	%
% (i.e. an intial condition file) for use by the DG_WASUPP pressure-Poisson	%
% model.									%
%										%
% The following is a list of the component that must be changed for each new	%
% hotstart file.								%
%										%
% hotrun : Name for hotstart run, will be written on the first line of the file %
% hottype : Type of hotstart file (NODAL or MODAL).  Currently only "Nodal" is	%
%           compatible so this shouldn't be changed.				%
% hotfile : Name of actual file, usually *****.67				%
% gridfile : Gridfile to be used for hotstart.					%
% h0 : If initial dry states are anticipated, this is the minimum water level.	%
% p : The order of the hotstart file, p>=1.  p=1 is standard for linear.	%
% zfun : A function of x that represents the free-surface displacement		%
% ufun : A function of x that represents the velocity				%
%										%
% There are a few pre-loaded hotstart setups which can be uncommented.		%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Hotstart variables
% hotrun = 'test';
% hottype = 'NODAL';
% hotfile = 'fort.67';
% h0 = 0.00001;
% p = 1;
% gridfile = 'fort.14';

% -- Solitary NonBreaking
%  hotrun = 'Solitary_Wave_NonBreaking';
%  hottype = 'NODAL';
%  hotfile = 'Solitary_nonbreak.67';
%  h0 = 0.00001;
%  p = 1;
%  gridfile = 'Solitary.14';
%  a = 0.0185;
%  hb = 0.4;
%  x0 = -15;
%  x1 = 1000;
%  c = sqrt(9.81*(hb+a));
%  kap = sqrt(3*a)/(2*hb*sqrt(hb+a));
%  zfun = @(x)a*sech(kap*(x-x0)).^2.*(x<=x1) - 10*(x>x1);
%  ufun = @(x)c*(1-hb./(zfun(x)+hb));

% -- Solitary Breaking
%  hotrun = 'Solitary_Wave_Breaking';
%  hottype = 'NODAL';
%  hotfile = 'Solitary_break.67';
%  h0 = 0.00001;
%  p = 1;
%  gridfile = 'Solitary.14';
%  a = 0.3;
%  hb = 0.4;
%  x0 = -15;
%  x1 = 1000;
%  c = sqrt(9.81*(hb+a));
%  kap = sqrt(3*a)/(2*hb*sqrt(hb+a));
%  zfun = @(x)a*sech(kap*(x-x0)).^2.*(x<=x1) - 10*(x>x1);
%  ufun = @(x)c*(1-hb./(zfun(x)+hb));

% -- Hsiao Case 1
 hotrun = 'Hsiao_Case_1';
 hottype = 'NODAL';
 hotfile = 'Hsiao_case1.67';
 h0 = 0.00001;
 p = 1;
 gridfile = 'Hsiao_case1.14';
 a = 0.07;
 hb = 0.2;
 x0 = 3;
 x1 = 13.9;
 c = sqrt(9.81*(hb+a));
 kap = sqrt(3*a)/(2*hb*sqrt(hb+a));
 zfun = @(x)a*sech(kap*(x-x0)).^2.*(x<=x1) - 10*(x>x1);
 ufun = @(x)c*(1-hb./(zfun(x)+hb));

% -- Hsiao Case 2
%  hotrun = 'Hsiao_Case_2';
%  hottype = 'NODAL';
%  hotfile = 'Hsiao_case2.67';
%  h0 = 0.00001;
%  p = 1;
%  gridfile = 'Hsiao_case2.14';
%  a = 0.0638;
%  hb = 0.22;
%  x0 = 3;
%  x1 = 13.9;
%  c = sqrt(9.81*(hb+a));
%  kap = sqrt(3*a)/(2*hb*sqrt(hb+a));
%  zfun = @(x)a*sech(kap*(x-x0)).^2.*(x<=x1) - 10*(x>x1);
%  ufun = @(x)c*(1-hb./(zfun(x)+hb));

% -- Hsiao Case 3
%  hotrun = 'Hsiao_Case_3';
%  hottype = 'NODAL';
%  hotfile = 'Hsiao_case3.67';
%  h0 = 0.00001;
%  p = 1;
%  gridfile = 'Hsiao_case3.14';
%  a = 0.0589;
%  hb = 0.18+0.076;
%  x0 = 3;
%  x1 = 13.9;
%  c = sqrt(9.81*(hb+a));
%  kap = sqrt(3*a)/(2*hb*sqrt(hb+a));
%  zfun = @(x)a*sech(kap*(x-x0)).^2.*(x<=x1) - 10*(x>x1);
%  ufun = @(x)c*(1-hb./(zfun(x)+hb));

% -- Carrier and Greenspan
%  hotrun = 'Carrier_Greenspan';
%  hottype = 'NODAL';
%  hotfile = 'carrier.67';
%  h0 = 0.00001;
%  p = 1;
%  gridfile = 'carrier.14';
%  a = 0.0185;
%  hb = 1;
%  x0 = -30;
%  x1 = 9999;
%  c = sqrt(9.81*(hb+a));
%  kap = sqrt(3*a)/(2*hb*sqrt(hb+a));
%  zfun = @(x)a*sech(kap*(x-x0)).^2.*(x<=x1) - 10*(x>x1);
%  ufun = @(x)c*(1-hb./(zfun(x)+hb));

% -- Dambreak
%  hotrun = 'Dambreak';
%  hottype = 'NODAL';
%  hotfile = 'dambreak.67';
%  h0 = 0.00001;
%  p = 1;
%  gridfile = 'dambreak.14';
%  zfun = @(x)10*(x<0)+0*(x>=0);
%  ufun = @(x)0*x;

% -- Riemann
%  hotrun = 'Riemann';
%  hottype = 'NODAL';
%  hotfile = 'riemann.67';
%  h0 = 0.00001;
%  p = 1;
%  gridfile = 'riemann.14';
%  zfun = @(x)5*(x<=0)+10*(x>0);
%  ufun = @(x)40*(x>0);

% -- Parabolic Bowl
%  hotrun = 'Parabolic';
%  hottype = 'NODAL';
%  hotfile = 'parabolic.67';
%  h0 = 0.00001;
%  p = 1;
%  gridfile = 'parabolic.14';
%  xlen = 8000;
%  zfun = @(x)0.001*cos(2*pi/xlen*x);
%  ufun = @(x)0*x;

%% Read grid file
fid = fopen(gridfile);
tmp = textscan(fid,'%f %f',1,'headerlines',1);
nn = tmp{2};
tmp = textscan(fid,'%f %f %f',nn);

xt = tmp{2};
dt = tmp{3};

if p == 0
    x = (xt(2:end)+xt(1:end-1))/2;
    d = (dt(2:end)+dt(1:end-1))/2;
end

if strcmp(hottype,'NODAL')
    nnhot = p*(nn-1)+1;
    x = zeros(nnhot,1);
    d = x;
    dx = 2/p;
    for l = 1:nn-1
        le = xt(l+1)-xt(l);
        for i = 1:p+1
            loc = (l-1)*p+i;
            x(loc) = xt(l)+le/2*( (i-1)*dx );
            d(loc) = dt(l)+(i-1)*dx/2*(dt(l+1)-dt(l));
        end
    end    
end

%% Hotstart z and q

ze = zfun(x);
% Make sure that the free-surface is beneath the bathymetry
for i = 1:length(ze)
    if (ze(i)+d(i))<=0
        ze(i) = h0-d(i);
    end
end 
qe = ufun(x).*(ze+d);

%% Write hotstart file
fid  = fopen(hotfile,'w');
fprintf(fid,'%s\n',hotrun);
fprintf(fid,'%s\n',hottype);
fprintf(fid,'%d\n',nnhot);
for l = 1:nnhot
    fprintf(fid,'%f %f %f\n',[x(l),ze(l),qe(l)]);
end
fclose(fid);

%% Visualize solution
figure(2)
plot(x,ze,xt,-dt,'--',x,qe,'m')
