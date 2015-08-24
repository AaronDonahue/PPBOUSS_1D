%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Matlab tool to visualize results					%
%										%
% This Matlab code has been designed to visualize results of a DG_WASUPP model  %
% run.										%
%										%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Determine which files to visualize
[file14, path14] = uigetfile('*.14','Pick the grid file');
[file63, path63] = uigetfile('*.63','Pick the global output file');
[file64, path64] = uigetfile('*.64','Pick the global output file');

%% load grid file
tmppath = sprintf('%s%s',path14,file14);
fid = fopen(tmppath);
tmp = textscan(fid,'%f %f',1,'headerlines',1);
ne = tmp{1};
nn = tmp{2};
tmp = textscan(fid,'%f %f %f',nn);
x = tmp{2};
d = tmp{3};
fclose(fid);

%% load output file
tmppath63 = sprintf('%s%s',path63,file63);
tmpeta = dlmread(tmppath63);
p = tmpeta(1,2);
if ne~=tmpeta(1,1)
    fprintf('Output file doesn''t match grid file!\n')
    return
end
tmppath64 = sprintf('%s%s',path64,file64);
tmpq = dlmread(tmppath63);
if p~=tmpq(1,2) || ne~=tmpq(1,1)
    fprintf('Output file doesn''t match 63 file!\n')
    return
end

t = tmpeta(2:ne+1:end,1);
eta = cell(length(t),1);
q = eta;
wd = eta;

%% DG Modes (NODAL)
if p == 1
    phi{1} = @(x)1/2*(1-x);
    phi{2} = @(x)1/2*(1+x);
elseif p == 2
    phi{1} = @(x)1/2*(x-1).*x;
    phi{2} = @(x)1-x.^2;
    phi{3} = @(x)1/2*(1+x).*x;
elseif p == 3
    phi{1} = @(x)1/16*(-1+x+9*x.^2-9*x.^3);
    phi{2} = @(x)9/16*(1-3*x-x.^2+3*x.^3);
    phi{3} = @(x)-9/16*(-1-3*x+x.^2+3*x.^3);
    phi{4} = @(x)1/16*(-1-x+9*x.^2+9*x.^3);
end

%% Read in data
pt = linspace(-1,1,p+1);
for i = 1:length(t)
    eta{i} = zeros(ne,p+1);
    q{i} = zeros(ne,p+1);
    wd{i} = zeros(ne,1);
    loc = (i-1)*(ne+1)+2;
    try
        eta{i}(:,:) = tmpeta(loc+1:loc+ne,2:2+p);
        wd{i}(:,:)  = tmpeta(loc+1:loc+ne,end);
    catch
        eta{i}(:,:) = 0;
    end
    try
        q{i}(:,:) = tmpq(loc+1:loc+ne,2:2+p);
    catch
        q{i}(:,:) = 0;
    end    
    % build solution from modes    
    etmp = eta{i};
    eta{i} = 0*eta{i};
    for j = 1:p+1
        for k = 1:p+1
            eta{i}(:,j) = eta{i}(:,j) + etmp(:,k)*phi{k}(pt(j));
        end
    end
    %
    eta{i} = eta{i}';
    
    qtmp = q{i};
    q{i} = 0*q{i};
    for j = 1:p+1
        for k = 1:p+1
            q{i}(:,j) = q{i}(:,j) + qtmp(:,k)*phi{k}(pt(j));
        end
    end
    %
    q{i} = q{i}';
    
end

xdg = zeros(ne,p+1);
for i = 1:ne
    xdg(i,:) = linspace(x(i),x(i+1),p+1);
end
xdg = xdg';

%% Plot solution
figure(1)
clf
for i = 1%:length(t)
    plot(xdg,eta{i},'b')
    hold on
    plot(xdg,q{i},'r--')
    plot(x,-d,'y')
    hold off
%     xlim([12 16])
    title(t(i))
    pause(0.001)    
end