clc
clear

% Set parameters
eps          = 5;
rbf_base_opt = 'PHS'; % Choose from PHS, GA
res          = pi/40;
res_plot     = pi/100;

% Calculate grids
x1d          = -pi:res:pi;
y1d          = -pi:res:pi;

nx = size(x1d,2);
ny = size(y1d,2);

[x2d,y2d] = meshgrid(x1d,y1d);

x = reshape(x2d,[],1);
y = reshape(y2d,[],1);

center_point = find(x==0&y==0);

f               = zeros(size(x));
f(center_point) = 1;

% Calculate weights matrix
r = distance(0,0,x,y,'radians');
K = rbf_base(r,eps,rbf_base_opt); % Calculate kernal function
A = bsxfun(@circshift,K,0:size(K,1)-1); % extent K to a 2d matrix
w = A \ f;

% Interpolate
x1d_plot = -pi:res_plot:pi;
y1d_plot = -pi:res_plot:pi;

[x2d_plot,y2d_plot] = meshgrid(x1d_plot,y1d_plot);

x_plot = reshape(x2d_plot,[],1);
y_plot = reshape(y2d_plot,[],1);

nx_plot = size(x1d_plot,2);
ny_plot = size(y1d_plot,2);
n_plot  = nx_plot*ny_plot;

f_plot_1d = zeros(n_plot,1);
for i = 1:n_plot
    r_plot       = distance(x_plot(i),y_plot(i),x,y,'radians');
    phi_plot     = rbf_base(r_plot,eps,rbf_base_opt);
    f_plot_1d(i) = sum(w .* phi_plot);
end

f_plot = reshape(f_plot_1d,nx_plot,ny_plot);

pic = surf(f_plot);
set(pic,'edgecolor','none')
colormap(jet)




%%
function phi = rbf_base(r,eps,opt)
if strcmp(opt,'GA')
    phi = exp(-(eps * r).^2);
elseif strcmp(opt,'PHS')
    phi = r.^eps;
end
end