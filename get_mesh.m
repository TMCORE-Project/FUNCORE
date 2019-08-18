function mesh = get_mesh(mesh_file,nSamples)

% Define Constants
mesh.omega = 7.292*10^-5;        %2.0 * pi / 86400.0; 
mesh.a     = 6371220.0;          %6371229.0
mesh.g     = 9.80616;
mesh.d2r   = pi/180;
mesh.r2d   = 180/pi;

xCell   = ncread(mesh_file,'xCell');
yCell   = ncread(mesh_file,'yCell');
zCell   = ncread(mesh_file,'zCell');

lonCell = ncread(mesh_file,'lonCell');
latCell = ncread(mesh_file,'latCell');

nCells  = size(xCell,1);

mesh.xCell   = xCell;
mesh.yCell   = yCell;
mesh.zCell   = zCell;
mesh.lonCell = lonCell;
mesh.latCell = latCell;
mesh.nCells  = nCells;

mesh.sinlon  = sin(mesh.lonCell);
mesh.coslon  = cos(mesh.lonCell);
mesh.sinlat  = sin(mesh.latCell);
mesh.coslat  = cos(mesh.latCell);
mesh.f       = 2 * mesh.omega * mesh.sinlat;

mesh.coord   = [xCell,yCell,zCell];

kdtree_searcher = KDTreeSearcher(mesh.coord);%,'BucketSize',100);
mesh.kdtree     = knnsearch(kdtree_searcher,mesh.coord,'K',nSamples);

distance = zeros(nCells,nSamples);
for iCell = 1:nCells
    chosen_coord     = mesh.coord(mesh.kdtree(iCell,:),:);
    origin_coord     = mesh.coord(mesh.kdtree(iCell,1),:);
    coord_diff       = chosen_coord - origin_coord;
    distance(iCell,:)= sqrt(sum(coord_diff.^2,2));
end

mesh.distance = distance;

mesh.weights = repmat(mesh.nCells/(4*pi),[1 mesh.nCells]); % Extract the quadrature weights.

% Variables for projecting an arbitrary Cartesian vector onto the surface
% of the sphere.
x2 = xCell.^2; xy = xCell.*yCell;
y2 = yCell.^2; xz = xCell.*zCell;
z2 = zCell.^2; yz = yCell.*zCell;

mesh.p_u = [1-x2  -xy   -xz];
mesh.p_v = [-xy  1-y2   -yz];
mesh.p_w = [-xz   -yz  1-z2];

% vectors for translating the field in Cartesian coordinates to a field
% in spherical coordinates.
mesh.c2s_u = [-sin(lonCell) -sin(latCell).*cos(lonCell)];
mesh.c2s_v = [cos(lonCell) -sin(latCell).*sin(lonCell)];
mesh.c2s_w = [zeros(size(lonCell)) cos(latCell)];

% Transformation for converting the latitudinal velocity to cartesian velocity.
mesh.s2c_u = [-sin(lonCell), cos(lonCell), zeros(size(lonCell))];
% Transformation for converting the logitudinal velocity to cartesian velocity.
mesh.s2c_v = [-cos(lonCell).*sin(latCell), -sin(lonCell).*sin(latCell), cos(latCell)];
