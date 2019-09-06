clc
clear

eps       = 9;
base_type = 2;

root_path = 'E:\Study\Models\FUNCORE';

% Choose mesh file
% mesh_file = 'x1.2562.grid.nc';
mesh_file = 'x1.40962.grid.nc';
% mesh_file = 'x1.163842.grid.nc';
% mesh_file = 'x1.655362.grid.nc';

nSamples = 400; % Number of sample points

d2r = pi/180;
r2d = 180/pi;

mesh_file = [root_path,'\',mesh_file];
mesh      = get_mesh(mesh_file,nSamples);

iCell    = 1;
dist     = distance(mesh.latCell(iCell),mesh.lonCell(iCell),mesh.latCell,mesh.lonCell,'radians');
[dist,I] = sort(dist,'ascend');
dist     = dist(1:nSamples);
I        = I   (1:nSamples);
lonCell  = mesh.lonCell(I);
latCell  = mesh.latCell(I);

r = zeros(nSamples,nSamples);
for i = 1:nSamples
    r(i,:) = distance(latCell(i),lonCell(i),latCell,lonCell,'radians');
end

% Select base function
if base_type == 1
    rbf_base = @(r,eps                    ) exp(-(eps.*r).^2);
    rbf_dlon = @(r,eps,lon,lat,lon_c,lat_c) -2*eps^2.*exp(-(eps.*r).^2).*cos(lat_c).*cos(lat).*sin(lon_c-lon);
    rbf_dlat = @(r,eps,lon,lat,lon_c,lat_c) -2*eps^2.*exp(-(eps.*r).^2).*(sin(lat_c).*cos(lat).*cos(lon_c-lon)-cos(lat_c).*cos(lat));
elseif base_type == 2
    rbf_base = @(r,eps                    ) r.^eps;
    rbf_dlon = @(r,eps,lon,lat,lon_c,lat_c) eps * r.^(eps-2).*cos(lat_c).*cos(lat).*sin(lon_c-lon);
    rbf_dlat = @(r,eps,lon,lat,lon_c,lat_c) eps * r.^(eps-2).*(sin(lat_c).*cos(lat).*cos(lon_c-lon)-cos(lat_c).*cos(lat));
elseif base_type == 3
    rbf_base = @(r,eps                    ) sqrt(1+(eps.*r).^2);
    rbf_dlon = @(r,eps,lon,lat,lon_c,lat_c) eps.^2 ./ sqrt(1+(eps*r).^2).*cos(lat_c).*cos(lat).*sin(lon_c-lon);
    rbf_dlat = @(r,eps,lon,lat,lon_c,lat_c) eps.^2 ./ sqrt(1+(eps*r).^2).*(sin(lat_c).*cos(lat).*cos(lon_c-lon)-cos(lat_c).*cos(lat));
end

K = rbf_base(r,eps);
% K = bsxfun(@circshift,K1d,0:size(K1d,1)-1); % extent K to a 2d matrix

[psi,beta] = orthogonalization(K);
beta       = beta';
orth_check = psi*psi';

dKdlon = rbf_dlon(r,eps,mesh.lonCell(I),mesh.latCell(I),mesh.lonCell(iCell),mesh.latCell(iCell));
dKdlat = rbf_dlat(r,eps,mesh.lonCell(I),mesh.latCell(I),mesh.lonCell(iCell),mesh.latCell(iCell));

Llon = dKdlon*beta'*beta;
Llat = dKdlat*beta'*beta;

p = pcolor(flipud(psi));
set(p,'edgecolor','none')
colormap(jet)
