clc
clear

% Choose mesh file
mesh_file = 'x1.2562.grid.nc';
% mesh_file = 'x1.40962.grid.nc';
% mesh_file = 'x1.163842.grid.nc';
% mesh_file = 'x1.655362.grid.nc';

% Define time(seconds)
run_day          = 1;
run_hour         = 0;
run_minute       = 0;
run_second       = 0;
time_step        = 600;
history_interval = 600;
temporal_scheme  = 'RK4';

% Select case
case_type       = 6;

nSample           = 31; % Number of sample points
base_opt          = 2; % 1 for r^m base, 2 for Gaussian base
viscosity_stencil = nSample;
viscosity_order   = 4;
viscosity_coef    = 1/6371229^40;
shape_param       = 20;  % shape parameter for RBF
poly_order        = 1;   % order of polynominal

mesh = get_mesh(mesh_file,nSample);

% Calculate weight matrix
w = zeros(mesh.nCells,nSample,3);
parfor iCell = 1: mesh.nCells
    disp(['Calculate weights for cell # ',num2str(iCell)])
    
    x   = mesh.xCell  (mesh.kdtree(iCell,:));
    y   = mesh.yCell  (mesh.kdtree(iCell,:));
    z   = mesh.zCell  (mesh.kdtree(iCell,:));
    
    w(iCell,:,:) = gen_weights(x,y,z,shape_param,poly_order,base_opt);
end

mesh.dx = squeeze(w(:,:,1));
mesh.dy = squeeze(w(:,:,2));
mesh.dz = squeeze(w(:,:,3));

% Calculate Laplacian operator
[~, ~, ~, ~, ~, mesh.L] = viscosity(mesh.coord,mesh.kdtree,shape_param,viscosity_stencil,viscosity_order,viscosity_coef,3,mesh.a);

% Initial field
[stat,mesh] = test_case(mesh,case_type);

run_time = run_day*86400 + run_hour*3600 + run_minute*60 + run_second;
run_step = run_time / time_step;

for it = 1:run_step
    disp(['step ',num2str(it),'/',num2str(run_step)])
    [stat] = temporal_integration(stat,mesh,time_step,temporal_scheme);
    
    integral_time = it*time_step;
    
    % Output
    if rem(integral_time,history_interval)==0 && integral_time>=history_interval
        % plot initial field for check
        lon1d = -180:1:180;
        lat1d = -90:1:90;
        
        [lon2d,lat2d] = meshgrid(lon1d,lat1d);
        
        var      = stat.gh;
        var_plot = griddata(mesh.lonCell*mesh.r2d,mesh.latCell*mesh.r2d,var,lon2d,lat2d);
        
        pcolor(lon2d,lat2d,var_plot);
        shading interp
        colormap(jet)
        colorbar
        print(gcf,'-r300','-dpng',['rbf',num2str(it,'%05d'),'.png']);
    end
end
