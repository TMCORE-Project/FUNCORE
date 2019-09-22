clc
clear

eps       = 7;

% Select base function
base_type = 2;

plot_res  = 0.5;

wave_number = 1;

root_path = 'E:\Study\Models\FUNCORE';

% Choose mesh file
mesh_file = 'x1.2562.grid.nc';
% mesh_file = 'x1.2562.grid.modified.nc';
% mesh_file = 'x1.10242.grid.nc';
% mesh_file = 'x1.40962.grid.nc';
% mesh_file = 'x1.163842.grid.nc';
% mesh_file = 'x1.655362.grid.nc';

nSamples = 2562; % Number of sample points

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
areaCell = mesh.areaCell;

r = zeros(nSamples,nSamples);
for i = 1:nSamples
    r(i,:) = distance(latCell(i),lonCell(i),latCell,lonCell,'radians');
end

% r = 2 * tan(r/2);

K = rbf_base(r,eps,base_type);

dKdlon = zeros(nSamples,nSamples);
dKdlat = zeros(nSamples,nSamples);
for i = 1:nSamples
    dKdlon(i,:) = rbf_dlon(r(i,:),eps,lonCell',latCell',lonCell(i),latCell(i),base_type);
    dKdlat(i,:) = rbf_dlat(r(i,:),eps,lonCell',latCell',lonCell(i),latCell(i),base_type);
end

Dlon = dKdlon / K;
Dlat = dKdlat / K;

f        = sin(wave_number*lonCell) + cos(wave_number*lonCell) + sin(wave_number*latCell) + cos(wave_number*latCell);
dfdlon_a = wave_number*cos(wave_number*lonCell) - wave_number*sin(wave_number*lonCell);
dfdlat_a = wave_number*cos(wave_number*latCell) - wave_number*sin(wave_number*latCell);
dfdlon_n = Dlon * f;
dfdlat_n = Dlat * f;

L2lon = sqrt(sum((dfdlon_n-dfdlon_a).^2.*areaCell) / sum(areaCell)) ...
      / sqrt(sum(dfdlon_a.^2.*areaCell) / sum(areaCell));

% Plot eigenvalue of differential operator
figure
lambda=eig(Dlon);
lambda_r = real(lambda);
lambda_i = imag(lambda);
plot(real(lambda),imag(lambda),'.')
spectrum_radius = max(sqrt(lambda_r.^2 + lambda_i.^2));

% Plot f
x = -pi:plot_res*d2r:pi;
y = -pi/2:plot_res*d2r:pi/2;

[x2d,y2d] = meshgrid(x,y);

f2d = griddata(lonCell,latCell,f,x2d,y2d);

figure
pcolor(x2d*r2d,y2d*r2d,f2d)
shading interp
colormap(jet)
colorbar

% Plot dfdlon, dfdlat analytical resolution
dfdlon2d = griddata(lonCell,latCell,dfdlon_a,x2d,y2d);
dfdlat2d = griddata(lonCell,latCell,dfdlat_a,x2d,y2d);

figure
pcolor(x2d*r2d,y2d*r2d,dfdlon2d);
set(gca,'clim',[-1.5 1.5])
shading interp
colormap(jet)
colorbar

figure
pcolor(x2d*r2d,y2d*r2d,dfdlat2d)
set(gca,'clim',[-1 1.5])
shading interp
colormap(jet)
colorbar

% Plot dfdlon, dfdlat numerical resolution
dfdlon2d = griddata(lonCell,latCell,dfdlon_n,x2d,y2d);
dfdlat2d = griddata(lonCell,latCell,dfdlat_n,x2d,y2d);

figure
pcolor(x2d*r2d,y2d*r2d,dfdlon2d)
set(gca,'clim',[-1.5 1.5])
shading interp
colormap(jet)
colorbar

figure
pcolor(x2d*r2d,y2d*r2d,dfdlat2d)
set(gca,'clim',[-1 1.5])
shading interp
colormap(jet)
colorbar

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

x      = cos(lon_c-lon).*cos(lat_c).*cos(lat)+sin(lat_c).*sin(lat);
xinv   = 1-x.^2;
xinv(abs(xinv)<10^-15)=0;
drdlon = sin(lon_c-lon).*cos(lat_c).*cos(lat)./xinv;

if base_type == 1
    dlon = -2*eps^2.*exp(-(eps.*r).^2).*r.*drdlon;
elseif base_type == 2
    dlon = eps * r.^(eps-1).*drdlon;
elseif base_type == 3
    dlon = eps.^2 ./ sqrt(1+(eps*r).^2).*r.*drdlon;
end

dlon(xinv==0)=0;
end

function dlat = rbf_dlat(r,eps,lon,lat,lon_c,lat_c,base_type)

x      = cos(lon_c-lon).*cos(lat_c).*cos(lat)+sin(lat_c).*sin(lat);
xinv   = 1-x.^2;
xinv(abs(xinv)<10^-15)=0;
drdlat = -(-cos(lon_c-lon).*sin(lat_c).*sin(lat)+cos(lat_c).*sin(lat))./xinv;

if base_type == 1
    dlat = -2*eps^2.*exp(-(eps.*r).^2).*r.*drdlat;
elseif base_type == 2
    dlat = eps * r.^(eps-1).*drdlat;
elseif base_type == 3
    dlat = eps.^2 ./ sqrt(1+(eps*r).^2).*r.*drdlat;
end

dlat(xinv==0)=0;
end