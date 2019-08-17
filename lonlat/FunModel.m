clc
clear

mesh_file = 'x1.40962.grid.modified.nc';
% mesh_file = 'x1.655362.grid.nc';
nSample    = 100; % Number of sample points

mesh = get_mesh(mesh_file,nSample);

shape_coef = 100/mesh.a;

w = zeros(mesh.nCells,nSample,2);
for iCell = 1: mesh.nCells
    x   = mesh.xCell  (mesh.kdtree(iCell,:));
    y   = mesh.yCell  (mesh.kdtree(iCell,:));
    z   = mesh.zCell  (mesh.kdtree(iCell,:));
    lon = mesh.lonCell(mesh.kdtree(iCell,:));
    lat = mesh.latCell(mesh.kdtree(iCell,:));
    
    weights = gen_weights(x,y,z,lon,lat,mesh.a,shape_coef,-1);
    w(iCell,:,:) = weights;
end

dlon = squeeze(w(:,:,1));
dlat = squeeze(w(:,:,2));