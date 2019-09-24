clc
clear

eps       = 3;

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
x        = mesh.xCell;
y        = mesh.yCell;
z        = mesh.zCell;
coord    = [x,y,z];
dist     = pdist(coord);
dist     = squareform(dist);
lonCell  = mesh.lonCell;
latCell  = mesh.latCell;

r = dist;

K = rbf_base(r,eps,base_type);
% K = bsxfun(@circshift,K1d,0:size(K1d,1)-1); % extent K to a 2d matrix

[psi,R]    = qr(K);
beta       = K/R/K;
% [psi,beta] = orthogonalization(K);
% beta       = beta';
orth_check = psi*psi';
orth_check(orth_check>0.99) = 0;

dKdlon = zeros(nSamples,nSamples);
dKdlat = zeros(nSamples,nSamples);
for i = 1:nSamples
    dKdlon(i,:) = rbf_dlon(r(i,:),eps,lonCell',latCell',lonCell(i),latCell(i),base_type);
    dKdlat(i,:) = rbf_dlat(r(i,:),eps,lonCell',latCell',lonCell(i),latCell(i),base_type);
end

Dlon = (beta'*beta*dKdlon)';
Dlat = (beta'*beta*dKdlat)';

f        = sin(wave_number*lonCell) + cos(wave_number*lonCell) + sin(wave_number*latCell) + cos(wave_number*latCell);
dfdlon_a = wave_number*cos(wave_number*lonCell) - wave_number*sin(wave_number*lonCell);
dfdlat_a = wave_number*cos(wave_number*latCell) - wave_number*sin(wave_number*latCell);
dfdlon_n = Dlon * f;
dfdlat_n = Dlat * f;

% Plot eigenvalue of differential operator
figure
lambda=eig(Dlon);
lambda_r = real(lambda);
lambda_i = imag(lambda);
plot(real(lambda),imag(lambda),'.')
spectrum_radius = max(sqrt(lambda_r.^2 + lambda_i.^2));

% p = pcolor(flipud(psi));
% set(p,'edgecolor','none')
% colormap(jet)

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

drdlat = sin(lat_c).*cos(lat).*cos(lon_c-lon) - cos(lat_c).*sin(lat);

if base_type == 1
    dlat = -2*eps^2.*exp(-(eps.*r).^2).*drdlat;
elseif base_type == 2
    dlat = eps * r.^(eps-2).*drdlat;
elseif base_type == 3
    dlat = eps.^2 ./ sqrt(1+(eps*r).^2).*drdlat;
end
end