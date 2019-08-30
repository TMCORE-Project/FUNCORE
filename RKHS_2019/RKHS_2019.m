clc
clear

eps       = 5;
base_type = 2;

root_path = 'E:\Study\Models\FUNCORE';

% Choose mesh file
mesh_file = 'x1.2562.grid.nc';
% mesh_file = 'x1.40962.grid.nc';
% mesh_file = 'x1.163842.grid.nc';
% mesh_file = 'x1.655362.grid.nc';

nSamples = 400; % Number of sample points

d2r = pi/180;
r2d = 180/pi;

mesh_file = [root_path,'\',mesh_file];
mesh      = get_mesh(mesh_file,nSamples);

iCell    = 100;
dist     = distance(mesh.lonCell(iCell),mesh.latCell(iCell),mesh.lonCell,mesh.latCell,'radians');
[dist,I] = sort(dist,'ascend');
dist     = dist(1:nSamples);
I        = I   (1:nSamples);

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

r = rbf_base(dist,eps);
K = r;
K = bsxfun(@circshift,K,0:size(K,1)-1); % extent K to a 2d matrix

[psi,beta] = orthogonalization(K);

dKdlon = rbf_dlon(K,eps,mesh.lonCell(I),mesh.latCell(I),mesh.lonCell(iCell),mesh.latCell(iCell));
dKdlat = rbf_dlat(K,eps,mesh.lonCell(I),mesh.latCell(I),mesh.lonCell(iCell),mesh.latCell(iCell));
[dpsidlon,dpsidlat] = Lpsi(beta,dKdlon,dKdlat);

Llon = dpsidlon*beta;
Llat = dpsidlat*beta;

% p = pcolor(flipud(psi));
% set(p,'edgecolor','none')
% colormap(jet)
