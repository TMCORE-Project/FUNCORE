function mesh = get_mesh(mesh_file,nSample)

% Define Constants
mesh.Omega = 7.292*10^-5;        %2.0 * pi / 86400.0; 
mesh.a     = 6371220.0;          %6371229.0
mesh.g     = 9.80616;

xCell   = ncread(mesh_file,'xCell');
yCell   = ncread(mesh_file,'yCell');
zCell   = ncread(mesh_file,'zCell');

lonCell = ncread(mesh_file,'lonCell');
latCell = ncread(mesh_file,'latCell');

nCells  = size(xCell,1);

mesh.xCell   = xCell * mesh.a;
mesh.yCell   = yCell * mesh.a;
mesh.zCell   = zCell * mesh.a;
mesh.lonCell = lonCell;
mesh.latCell = latCell;
mesh.nCells  = nCells;

mesh.sinlat  = sin(mesh.latCell);
mesh.f       = 2 * mesh.Omega * mesh.sinlat;

mesh.coord   = [xCell,yCell,zCell];

kdtree_searcher = KDTreeSearcher(mesh.coord);%,'BucketSize',100);
mesh.kdtree     = knnsearch(kdtree_searcher,mesh.coord,'K',nSample);