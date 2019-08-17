function mesh = get_mesh(mesh_file,nSamples)

% Define Constants
mesh.Omega = 7.292*10^-5;        %2.0 * pi / 86400.0; 
mesh.a     = 6371220.0;          %6371229.0
mesh.g     = 9.80616;

xCell   = ncread(mesh_file,'xCell') * mesh.a;
yCell   = ncread(mesh_file,'yCell') * mesh.a;
zCell   = ncread(mesh_file,'zCell') * mesh.a;

lonCell = ncread(mesh_file,'lonCell');
latCell = ncread(mesh_file,'latCell');

nCells  = size(xCell,1);

mesh.xCell   = xCell;
mesh.yCell   = yCell;
mesh.zCell   = zCell;
mesh.lonCell = lonCell;
mesh.latCell = latCell;
mesh.nCells  = nCells;

mesh.sinlat  = sin(mesh.latCell);
mesh.f       = 2 * mesh.Omega * mesh.sinlat;

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