clc
clear

% Choose mesh file
mesh_file = 'x1.2562.grid.nc';
% mesh_file = 'x1.40962.grid.nc';
% mesh_file = 'x1.163842.grid.nc';
% mesh_file = 'x1.655362.grid.nc';

% Choose output path
output_path = '.\picture';

% Define time(seconds)
run_day          = 33;
run_hour         = 0;
run_minute       = 0;
run_second       = 0;
time_step        = 600;
history_interval = 3600;
temporal_scheme  = 'RK4';

% Select case
case_type         = 6;

nSample           = 101; % Number of sample points
base_opt          = 2; % 1 for r^m base, 2 for Gaussian base
viscosity_stencil = nSample;
viscosity_order   = 4;
viscosity_coef    = 1/6371229^40;
eps               = 7;  % shape parameter for RBF
poly_order        = 1;   % order of polynominal

mesh = get_mesh(mesh_file,nSample);

% Calculate weight matrix
w = zeros(mesh.nCells,nSample,3);
parfor iCell = 1: mesh.nCells
% for iCell = 1: mesh.nCells
    disp(['Calculate weights for cell # ',num2str(iCell)])
    
    x   = mesh.xCell  (mesh.kdtree(iCell,:));
    y   = mesh.yCell  (mesh.kdtree(iCell,:));
    z   = mesh.zCell  (mesh.kdtree(iCell,:));
    
    xi = [x,y,z];
    xc = mesh.coord(iCell,:);
    
    w(iCell,:,:) = gen_weights(x,y,z,eps,poly_order,base_opt);
    dx(iCell,:) = rbfga_weights('x',eps,xi,xc);
    dy(iCell,:) = rbfga_weights('y',eps,xi,xc);
    dz(iCell,:) = rbfga_weights('z',eps,xi,xc);
end

mesh.dx = squeeze(dx(:,:));
mesh.dy = squeeze(dy(:,:));
mesh.dz = squeeze(dz(:,:));

% Calculate Laplacian operator
[~, ~, ~, ~, ~, mesh.L] = viscosity(mesh.coord,mesh.kdtree,eps,viscosity_stencil,viscosity_order,viscosity_coef,3,mesh.a);

% Initial field
[stat,mesh] = test_case(mesh,case_type);

run_time = run_day*86400 + run_hour*3600 + run_minute*60 + run_second;
run_step = run_time / time_step;

output_idx = 0;

for it = 1:run_step
    disp(['step ',num2str(it),'/',num2str(run_step)])
    [stat] = temporal_integration(stat,mesh,time_step,temporal_scheme);
    
    integral_time = it*time_step;
    
    % Output
    if rem(integral_time,history_interval)==0 && integral_time>=history_interval
        output_idx = output_idx + 1;
        
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
        title(['RBF model at hour',num2str(output_idx)])
        print(gcf,'-r300','-dpng',[output_path,'\','rbf',num2str(output_idx,'%05d'),'.png']);
    end
end
