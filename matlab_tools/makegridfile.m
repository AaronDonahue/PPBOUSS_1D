%% Matlab script to construct a grid file

runname = 'dambreak';
gridfile = 'dambreak.14';

% runname = 'riemann';
% gridfile = 'riemann.14';

% runname = 'parabolic';
% gridfile = 'parabolic.14';
%% Determine x vector
% nn = 501;
% x = linspace(-10,10,nn);
dx = 10;
x = -300:dx:300;
% dx = 1.25;
% x = -200:dx:400;
% dx = 12.5;
% x = -4000:dx:4000;
nn = length(x);

% % make unstructured:
% for i = 2:length(x)-1
%     x(i) = chop(x(i) + (rand(1)-0.5)/4*dx,6);
% end

%% Determine depth
% dep_fun = @(x)-1.6*10^(-7)*x.^2;
dep_fun = @(x)0*x;
% dep_fun = @(x)0*x+1;
% dep_fun = @(x)0*x+1-2/10*x.*(x>0);
% dep_fun = @(x)0*x+1-(1+tanh((x-5)/2))/2/10;
dep = dep_fun(x);

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
