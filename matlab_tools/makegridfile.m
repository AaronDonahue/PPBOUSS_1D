%% Matlab script to construct a grid file

runname = 'Ting_Kirby';
gridfile = 'TingKirby.14';

% runname = 'SolitaryWave_break';
% gridfile = 'Solitary.14';

% runname = 'Hsiao_case3';
% gridfile = 'Hsiao_case3.14';

% runname = 'Carrier_Greenspan';
% gridfile = 'carrier.14';

% runname = 'dambreak';
% gridfile = 'dambreak.14';

% runname = 'riemann';
% gridfile = 'riemann.14';

% runname = 'parabolic';
% gridfile = 'parabolic.14';
%% Determine x vector
dx = 0.025;    % Solitary
x = -25:dx:25;
% dx = 0.02;    % Hsiao
% x = -10:dx:25;
% dx = 0.1;  % Carrier-Greenspan
% x = -45:dx:5;
% dx = 10;
% x = -300:dx:300;
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
% dep_fun = @(x)0*x;
% dep_fun = @(x)0*x+.4;
% dep_fun = @(x)1 - 1/19.85*(x+19.85).*(x>=-19.85); % Carrier Greenspan


% hb = 0.18;  % Hsiao
% eadj = 0.02;
% eadj = 0.04;
% eadj = 0.076;
% dep_fun = @(x) (hb - 1/20*(x-10).*(x>=10).*(x<=13.6) ...
%        - (1/4*(x-13.6)+hb).*(x>13.6).*(x<13.9) ...
%        - (0.076+hb).*(x>=13.9).*(x<=13.948) ...
%        - (-1/1.8*(x-13.948)+0.076+hb).*(x>=13.948).*(x<=14.045) ...
%        - 1/20*(x-10).*(x>14.045))+eadj;

% dep_fun = @(x)0*x+1-2/10*x.*(x>0);
% dep_fun = @(x)0*x+1-(1+tanh((x-5)/2))/2/10;

% hb = 0.3;  % Solitary Wave
% hb_0 = -hb*19.85;
% dep_fun = @(x)hb - 1/19.85*(x-hb_0).*(x>=hb_0);

hb = 0.4;
hb_0 = 0;
hb_slope = 1/35;
dep_fun = @(x)0*x+hb - .02/.2*(x-hb_0).*(x>=hb_0).*(x<hb_0+0.2)-(0.02+hb_slope*(x-hb_0-0.2)).*(x>=hb_0+0.2);

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
