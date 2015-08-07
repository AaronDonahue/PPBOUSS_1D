%% Script to construct hotstart file
hotrun = 'test';
hottype = 'NODAL';
hotfile = 'Solitary_nonbreak.67';
% hotfile = 'Solitary_break.67';
% hotfile = 'Hsiao_case1.67';
% hotfile = 'Hsiao_case2.67';
% hotfile = 'Hsiao_case3.67';
% hotfile = 'carrier.67';
% hotfile = 'fort.67';
% hotfile = 'dambreak.67';
% hotfile = 'riemann.67';
% hotfile = 'parabolic.67';

h0 = 0.00001;
p = 1;

if p == 0
    hottype = 'MODAL';
end
%% Read grid file
% gridfile = 'fort.14';
gridfile = 'Solitary.14';
% gridfile = 'Hsiao_case1.14';
% gridfile = 'Hsiao_case2.14';
% gridfile = 'Hsiao_case3.14';
% gridfile = 'carrier.14';
% gridfile = 'dambreak.14';
% gridfile = 'riemann.14';
% gridfile = 'parabolic.14';
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
% xlen = x(end)-x(1);
% zfun = @(x)0.001*cos(2*pi/xlen*x);
% zfun = @(x)0*x-d+10*(x<0)+0*(x>=0);
% zfun = @(x)0.1*exp(-x.^2);
% zfun = @(x)0*x+0.1;
% zfun = @(x)5*(x<=0)+10*(x>0);
% zfun = @(x)1/(1-0.41884)+1.6*10^(-7)*(0.41884^2-1)*x.^2/(1+0.41884)^2;
% ufun = @(x)0*x;
% ufun = @(x)40*(x>0);
 
 
% a = 0.0185; % Carrier_Greenspan
% hb = 1;
% x0 = -30;
% x1 = 9999;

% a = 0.07;  % Hsiao
% hb = 0.2;
% a = 0.0638;
% hb = 0.22;
% a = 0.0589;
% hb = 0.18+0.076;
% x0 = 3;
% x1 = 13.9;

a = 0.0185;
hb = 0.3;
x0 = -15;
x1 = 1000;

c = sqrt(9.81*(hb+a));
kap = sqrt(3*a)/(2*hb*sqrt(hb+a));

zfun = @(x)a*sech(kap*(x-x0)).^2.*(x<=x1) - 10*(x>x1);
ufun = @(x)c*(1-hb./(zfun(x)+hb));


ze = zfun(x);
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
