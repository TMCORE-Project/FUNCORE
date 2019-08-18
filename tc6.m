%
% Set up for the Williamson test case 6.
%
function [stat,mesh] = tc6(mesh)

% Extract out the constants from the atm structure to make the code easier
% to read. 
x = mesh.xCell;
y = mesh.yCell;
z = mesh.zCell;

a      = mesh.a;
omega  = mesh.omega;
omg    = 7.848e-6;  % A constant (1/s)
g      = mesh.g;    % Gravitational constant (m/s^2).
gh0    = g*8e3;     % Initial condition for the geopotential field (m^2/s^2).
K      = 7.848e-6;  % A constant (1/s)
R      = 4;         % Wave number of the R-H wave

lonCell = mesh.lonCell;
latCell = mesh.latCell;

us = [a*omg*cos(latCell)+a*K*cos(latCell).^(R-1).*(R*sin(latCell).^2-cos(latCell).^2).*cos(R*lonCell) ...
     -a*K*R*cos(latCell).^(R-1).*sin(latCell).*sin(R*lonCell)];

% Transformation for converting the latitudinal velocity to cartesian velocity.
s2c_u = mesh.s2c_u;
% Transformation for converting the logitudinal velocity to cartesian velocity.
s2c_v = mesh.s2c_v;
uc = s2c_u.*repmat(us(:,1),1,3) + s2c_v.*repmat(us(:,2),1,3);

A = 0.5*omg*(2*omega+omg)*cos(latCell).^2 + 0.25*K^2*cos(latCell).^(2*R).*...
    ((R+1)*cos(latCell).^2 + (2*R^2-R-2) - 2*R^2*cos(latCell).^(-2));
B = (2*(omega+omg)*K)*cos(latCell).^R.*((R^2+2*R+2)-(R+1)^2*cos(latCell).^2)/((R+1)*(R+2));
C = 0.25*K^2*cos(latCell).^(2*R).*((R+1)*cos(latCell).^2 - (R+2));

gh   = gh0 + a^2*A + a^2*B.*cos(R*lonCell) + a^2*C.*cos(2*R*lonCell); % height field without mean offset
ghs  = zeros(size(gh));
dghs = zeros(size(gh,1),3);

stat.gh   = gh;
stat.u    = uc(:,1);
stat.v    = uc(:,2);
stat.w    = uc(:,3);
stat.us   = us(:,1);
stat.vs   = us(:,2);
mesh.ghs  = ghs;
mesh.dghs = dghs;