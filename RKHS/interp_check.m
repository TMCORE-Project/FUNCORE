clc
clear

interp_resolution = 0.1;

%bf_opt --- base function option, 1 for r^eps, 2 for gaussian function exp(-(eps*r).^2)
bf_opt = 4;

% power of distance weight
eps = -4;

% interp domain scheme, choose from : 'fr' = fixed radius ;  'ng' = neighbours
r1 = 'fr';

% interp parameter, radius lenght if r1 == 'fr' / number of neighbours if  r1 =='ng'
r2 = 40;

x_dest = -10:interp_resolution:10;
y_dest = -10:interp_resolution:10;

x_src = -10:1:10;
y_src = -10:1:10;

nx = size(x_src,2);
ny = size(y_src,2);

[x2d,y2d] = meshgrid(x_src,y_src);

f0         = zeros(size(x2d));

f0((nx+1)/2,(ny+1)/2) = 1;

f = IDW(x2d,y2d,f0,x_dest,y_dest,eps,r1,r2,bf_opt);

% fig = pcolor(f);
% set(fig,'edgecolor','none')
% colormap(jet)

fig = surf(x_dest,y_dest,f);
% set(fig,'edgecolor','flat')
% shading interp
colormap(jet)