%
% Set up for the Williamson test case 5.
%
function [stat,mesh] = tc5(mesh)
a     = mesh.a;
omega = mesh.omega;
gh0   = mesh.g*5960;     % Initial condition for the geopotential field (m^2/s^2).
alpha = 0;           % Angle of rotation measured from the equator.
u0    = 20;          % Speed of rotation in meters/second

% Parameters for the mountain:
lambda_c = -pi/2;
theta_c  = pi/6;
mR       = pi/9;
hs0      = 2000;

% Compute the profile of the mountain (multiplied by gravity)
r2      = (mesh.lonCell-lambda_c).^2 + (mesh.latCell-theta_c).^2;
id      = r2 < mR^2;
ghs     = zeros(mesh.nCells,1);
ghs(id) = mesh.g*hs0*(1-sqrt(r2(id))/mR);
%%% Option: Gaussian mountain
% disp('Using Gaussian mountain.');
% % % atm.ghm = atm.g*atm.hm0*(exp(-2.25^2*(r2/atm.mR^2)));
% atm.ghm = atm.g*atm.hm0*(exp(-2.8*(r2/atm.mR^2)));

% Extract out the constants from the atm structure to make the code easier
% to read. 
x = mesh.xCell;
y = mesh.yCell;
z = mesh.zCell;

sinlat = mesh.sinlat;
coslat = mesh.coslat;
sinlon = mesh.sinlon;
coslon = mesh.coslon;

sin_alpha = sin(alpha);
cos_alpha = cos(alpha);

gh = gh0 - (a .* omega .* u0 + 0.5 .* u0^2) .* ( - coslon .* coslat .* sin_alpha + sinlat .* cos_alpha ).^2 - ghs;
uc = u0.*[-y.*cos_alpha, (x.*cos_alpha + z.*sin_alpha), -y.*sin_alpha];

% vectors for translating the field in Cartesian coordinates to a field
% in spherical coordinates.
c2s_u = mesh.c2s_u;
c2s_v = mesh.c2s_v;
c2s_w = mesh.c2s_w;

us(:,1) = c2s_u(:,1).*uc(:,1) + c2s_v(:,1).*uc(:,2) + c2s_w(:,1).*uc(:,3);
us(:,2) = c2s_u(:,2).*uc(:,1) + c2s_v(:,2).*uc(:,2) + c2s_w(:,2).*uc(:,3);

stat.gh   = gh;
stat.u    = uc(:,1);
stat.v    = uc(:,2);
stat.w    = uc(:,3);
stat.us   = us(:,1);
stat.vs   = us(:,2);
mesh.ghs  = ghs;

mesh.dghs(:,1) = deriv(mesh.ghs, mesh.dx, mesh.kdtree,a);
mesh.dghs(:,2) = deriv(mesh.ghs, mesh.dy, mesh.kdtree,a);
mesh.dghs(:,3) = deriv(mesh.ghs, mesh.dz, mesh.kdtree,a);
