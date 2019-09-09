clc
clear

eps       = 15;

% Select base function
base_type = 1;

root_path = 'E:\Study\Models\FUNCORE';

% Choose mesh file
mesh_file = 'x1.2562.grid.nc';
% mesh_file = 'x1.10242.grid.nc';
% mesh_file = 'x1.40962.grid.nc';
% mesh_file = 'x1.163842.grid.nc';
% mesh_file = 'x1.655362.grid.nc';

nSamples = 2562; % Number of sample points

d2r = pi/180;
r2d = 180/pi;

mesh_file = [root_path,'\',mesh_file];
mesh      = get_mesh(mesh_file,nSamples);

iCell    = 1000;
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

K = rbf_base(r,eps,base_type);

dKdlon = zeros(nSamples,nSamples);
dKdlat = zeros(nSamples,nSamples);
parfor i = 1:nSamples
    dKdlon(i,:) = rbf_dlon(r(i,:),eps,lonCell',latCell',lonCell(i),latCell(i),base_type);
    dKdlat(i,:) = rbf_dlat(r(i,:),eps,lonCell',latCell',lonCell(i),latCell(i),base_type);
end

Dlon = K \ dKdlon;
Dlat = K \ dKdlat;

f        = 1 + sin(lonCell) + cos(lonCell) + sin(latCell);
DfDlon_a = cos(lonCell) - sin(lonCell);
DfDlat_a = cos(latCell);
DfDlon_n = Dlon * f;
DfDlat_n = Dlat * f;

figure
plot(DfDlon_a,'b')
hold on
plot(DfDlon_n,'r')

% figure
% plot(DfDlat_a,'b')
% hold on
% plot(DfDlat_n,'r')

function rbf = rbf_base(r,eps,base_type)
if base_type == 1
    rbf = exp(-(eps.*r).^2);
elseif base_type == 2
    rbf = r.^eps;
elseif base_type == 3
    rbf = sqrt(1+(eps.*r).^2);
end
end

function dlon = rbf_dlon(r,eps,lon,lat,lon_c,lat_c,base_type)

drdlon = cos(lat_c).*cos(lat).*sin(lon_c-lon);

if base_type == 1
    dlon = -2*eps^2.*exp(-(eps.*r).^2).*drdlon;
elseif base_type == 2
    dlon = eps * r.^(eps-2).*drdlon;
elseif base_type == 3
    dlon = eps.^2 ./ sqrt(1+(eps*r).^2).*drdlon;
end
end

function dlat = rbf_dlat(r,eps,lon,lat,lon_c,lat_c,base_type)

drdlat = sin(lat_c).*cos(lat).*cos(lon_c-lon) - cos(lat_c).*cos(lat);

if base_type == 1
    dlat = -2*eps^2.*exp(-(eps.*r).^2).*drdlat;
elseif base_type == 2
    dlat = eps * r.^(eps-2).*drdlat;
elseif base_type == 3
    dlat = eps.^2 ./ sqrt(1+(eps*r).^2).*drdlat;
end
end