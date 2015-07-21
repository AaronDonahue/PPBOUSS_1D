%% Matlab tool to visualize results

[file14, path14] = uigetfile('*','Pick the grid file');
[file63, path63] = uigetfile('*','Pick the global output file');

% load grid file
tmppath = sprintf('%s%s',path14,file14);
fid = fopen(tmppath);
tmp = textscan(fid,'%f %f',1,'headerlines',1);
ne = tmp{1};
nn = tmp{2};
tmp = textscan(fid,'%f %f %f',nn);
x = tmp{2};
d = tmp{3};
fclose(fid);

% load output file
tmppath = sprintf('%s%s',path63,file63);
tmp = dlmread(tmppath);
if nn~=tmp(1,1)
    fprintf('Output file doesn''t match grid file!\n')
    return
end
t = tmp(2:nn+1:end,1);
eta = zeros(nn,length(t));
for i = 1:length(t)
    loc = (i-1)*(nn+1)+2;
    eta(:,i) = tmp(loc+1:loc+nn,2);
end

%% Plot solution
figure(1)
clf
for i = 1:length(t)
    plot(x,eta(:,i))
    hold on
    plot(x,-d,'y')
    hold off
%     xlim([x(1) x(end)]);
%     ylim([min(eta(:)) max(eta(:))])
%     xlim([4 6])
%     ylim([-0.02 0.2])
    title(t(i))
    pause(0.01)
end