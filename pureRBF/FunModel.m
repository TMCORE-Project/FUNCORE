clc
clear

mesh_file = '../x1.2562.grid.nc';
% mesh_file = 'x1.655362.grid.nc';
nSample    = 300; % Number of sample points

mesh = get_mesh(mesh_file,nSample);

shape_coef = 20;

w = zeros(mesh.nCells,nSample,3);
parfor iCell = 1: mesh.nCells
    disp(['Calculate weights for cell # ',num2str(iCell)])
    
    x   = mesh.xCell  (mesh.kdtree(iCell,:));
    y   = mesh.yCell  (mesh.kdtree(iCell,:));
    z   = mesh.zCell  (mesh.kdtree(iCell,:));
    lon = mesh.lonCell(mesh.kdtree(iCell,:));
    lat = mesh.latCell(mesh.kdtree(iCell,:));
    
    w(iCell,:,:) = gen_weights(x,y,z,shape_coef);
end

dx = squeeze(w(:,:,1));
dy = squeeze(w(:,:,2));
dz = squeeze(w(:,:,3));