clc
clear

% Set parameters
eps          = 5;
rbf_base_opt = 'PHS'; % Choose from PHS, GA
res          = 1/25;
res_plot     = 1/400;
x_start      = -1;
x_end        = 1;
y_start      = -1;
y_end        = 1;

% Calculate grids
x1d          = x_start:res:x_end;
y1d          = y_start:res:y_end;

nx = size(x1d,2);
ny = size(y1d,2);
n  = nx*ny;

[x2d,y2d] = meshgrid(x1d,y1d);

x     = reshape(x2d,[],1);
y     = reshape(y2d,[],1);
coord = [x,y];

center_point = find(x==0&y==0);

f               = zeros(size(x));
f(center_point) = 1;

% Calculate weights matrix
r = pdist(coord);
K = rbf_base(r,eps,rbf_base_opt); % Calculate kernal function
A = squareform(K);
w = A \ f;

% Interpolate
x1d_plot = x_start:res_plot:x_end;
y1d_plot = y_start:res_plot:y_end;

[x2d_plot,y2d_plot] = meshgrid(x1d_plot,y1d_plot);

x_plot = reshape(x2d_plot,[],1);
y_plot = reshape(y2d_plot,[],1);

nx_plot = size(x1d_plot,2);
ny_plot = size(y1d_plot,2);
n_plot  = nx_plot*ny_plot;

f_plot_1d = zeros(n_plot,1);
parfor i = 1:n_plot
    r_plot       = sqrt( (x_plot(i) - x).^2 + (y_plot(i) - y).^2 );
    phi_plot     = rbf_base(r_plot,eps,rbf_base_opt);
    f_plot_1d(i) = sum(w .* phi_plot);
end

f_plot = reshape(f_plot_1d,nx_plot,ny_plot);

pic = surf(x2d_plot,y2d_plot,f_plot);
% set(pic,'edgecolor','none')
colormap(jet)




%%
function phi = rbf_base(r,eps,opt)
if strcmp(opt,'GA')
    phi = exp(-(eps * r).^2);
elseif strcmp(opt,'PHS')
    phi = r.^eps;
end
end