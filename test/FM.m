clc
clear

mesh_file = '../x1.2562.grid.nc';

radius = 1;
eps    = 8;

xCell = ncread(mesh_file,'xCell');
yCell = ncread(mesh_file,'yCell');
zCell = ncread(mesh_file,'zCell');

lonCell = ncread(mesh_file,'lonCell');
latCell = ncread(mesh_file,'latCell');

nCell = size(xCell,1);

coord(:,1) = xCell;
coord(:,2) = yCell;
coord(:,3) = zCell;

lonlat(:,1) = lonCell;
lonlat(:,2) = latCell;

r_1d   = pdist(coord);
r      = squareform(r_1d);
% phi    = (eps * r).^5;
% dphidr = 5 * eps^5 * r.^3; % dphidr = dphidr / r, to aviod divide 0
phi    = exp( - eps^2 * r.^2 );
dphidr = -2 * eps^2 .* phi; % dphidr = dphidr / r, to aviod divide 0

inv_phi = inv(phi);

drdlambda = zeros(nCell,nCell);
drdtheta  = zeros(nCell,nCell);
for i = 1:nCell
    for j = 1:nCell
        drdlambda(i,j) = cos(latCell(j)) .* cos(latCell(i)) .* sin(lonCell(j) - lonCell(i)); % divide r is removed to avoid divide 0;
        drdtheta (i,j) = ( sin(latCell(j)) * cos(latCell(i)) * cos(lonCell(j)-lonCell(i)) - cos(latCell(j)) * sin(latCell(i)) ); % divide r is removed to avoid divide 0;
    end
end

dphidlambda = dphidr .* drdlambda;
dphidtheta  = dphidr .* drdtheta;

ddlambda = dphidlambda / phi;
ddtheta  = dphidtheta  / phi;