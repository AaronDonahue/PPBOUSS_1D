%% Make nodal attributes file
nodefile = 'fort.13';
noderun = 'TingKirby';
p = 1;

%% Read grid file
gridfile = 'TingKirby.14';
% gridfile = 'Solitary.14';
% gridfile = 'Hsiao_case1.14';
% gridfile = 'carrier.14';
% gridfile = 'fort.14';
% gridfile = 'dambreak.14';
% gridfile = 'riemann.14';
% gridfile = 'parabolic.14';
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
attr(ind).on        = 0;
attr(ind).default   = 0.004;
attr(ind).numnondef = 0;
attr(ind).nondefvalue = [];
attr(ind).nondefelem  = [];

% Sponge Generation
ind = ind+1;
attr(ind).name    = 'sponge_generation_layer';
attr(ind).unit    = 'unitless';
attr(ind).on      = 1;
attr(ind).default = 0;
attr(ind).numnondef = 0;
attr(ind).nondefvalue = [];
attr(ind).nondefelem  = [];

% Sponge Absorbing
ind = ind+1;
attr(ind).name    = 'sponge_absorbing_layer';
attr(ind).unit    = 'unitless';
attr(ind).on      = 0;
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