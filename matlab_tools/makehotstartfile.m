%% Script to construct hotstart file
hotrun = 'test';
hottype = 'NODAL';
% hotfile = 'dambreak.67';
hotfile = 'riemann.67';
% hotfile = 'parabolic.67';

h0 = 0.00001;
p = 1;

if p == 0
    hottype = 'MODAL';
end
%% Read grid file
% gridfile = 'fort.14';
% gridfile = 'dambreak.14';
gridfile = 'riemann.14';
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
% zfun = @(x)0*x-d+10*(x<0)+0*(x>=0);
% zfun = @(x)0.05*exp(-x.^2);
% zfun = @(x)0*x+0.1;
zfun = @(x)5*(x<=0)+10*(x>0);
% zfun = @(x)1/(1-0.41884)+1.6*10^(-7)*(0.41884^2-1)*x.^2/(1+0.41884)^2;
% ufun = @(x)0*x;
ufun = @(x)40*(x>0);

ze = zfun(x)+d;
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
plot(x,ze,xt,-dt)