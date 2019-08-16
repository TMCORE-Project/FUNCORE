clc
clear

mesh_file = 'x1.2562.grid.nc';

radius = 1;
eps    = 10;

xCell = ncread(mesh_file,'xCell');
yCell = ncread(mesh_file,'yCell');
zCell = ncread(mesh_file,'zCell');

lonCell = ncread(mesh_file,'lonCell');
latCell = ncread(mesh_file,'latCell');

nCell = size(xCell,1);

coord(:,1) = xCell;
coord(:,2) = yCell;
coord(:,3) = zCell;

r_1d   = pdist(coord);
r      = squareform(r_1d);
% phi    = (eps * r).^5;
% dphidr = 5 * eps^5 * r.^4;
% phi    = exp( - eps^2 * r.^2 );
% dphidr = -2 * eps^2 * r .* phi;
phi    = exp( - eps^2 * r.^2 );
dphidr = -2 * eps^2 * r .* phi;

% inv_phi = inv(phi);

drdx      = zeros(nCell,nCell);
drdy      = zeros(nCell,nCell);
drdz      = zeros(nCell,nCell);
for i = 1:nCell
    for j = 1:nCell
        if i~=j
            drdx(i,j) = (xCell(j)-xCell(i)) / r(i,j);
            drdy(i,j) = (yCell(j)-yCell(i)) / r(i,j);
            drdz(i,j) = (zCell(j)-zCell(i)) / r(i,j);
        else
            drdx(i,j) = 1;
            drdy(i,j) = 1;
            drdz(i,j) = 1;
        end
    end
end

dphidx  = dphidr .* drdx;
dphidy  = dphidr .* drdy;
dphidz  = dphidr .* drdz;

ddx  = dphidx / phi;
ddy  = dphidy / phi;
ddz  = dphidz / phi;