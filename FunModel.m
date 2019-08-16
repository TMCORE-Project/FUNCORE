clc
clear

mesh_file = 'x1.40962.grid.modified.nc';
nSample   = 60; % Number of sample points

mesh = get_mesh(mesh_file,nSample);

for iCell = 1: mesh.nCells
    x   = mesh.xCell  (mesh.kdtree(iCell,:));
    y   = mesh.yCell  (mesh.kdtree(iCell,:));
    z   = mesh.zCell  (mesh.kdtree(iCell,:));
    lon = mesh.lonCell(mesh.kdtree(iCell,:));
    lat = mesh.latCell(mesh.kdtree(iCell,:));
    w = gen_weights(x,y,z,lon,lat,mesh.a,7,7);
end