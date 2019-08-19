clc
clear

d2r = pi/180;
r2d = 180/pi;

res1 = 1;
res1 = res1*d2r;
lon1 = 0:res1:2*pi-res1;
lon1 = lon1';
lat1 = zeros(size(lon1));

dist1 = distance(0.,0.,lat1,lon1,'radians');

res2 = 2;
res2 = res2*d2r;
lon2 = 0:res2:2*pi-res2;
lon2 = lon2';
lat2 = zeros(size(lon2));

dist2 = distance(0.,0.,lat2,lon2,'radians');

res3 = 3;
res3 = res3*d2r;
lon3 = 0:res3:2*pi-res3;
lon3 = lon3';
lat3 = zeros(size(lon3));

dist3 = distance(0.,0.,lat3,lon3,'radians');

beta1 = ncread('E:\Study\Models\SplineDyn\genMatrix\run\matrix_1p0_reduce1p0.nc','lon_matrix_beta');
beta2 = ncread('E:\Study\Models\SplineDyn\genMatrix\run\matrix_2p0_reduce1p0.nc','lon_matrix_beta');
beta3 = ncread('E:\Study\Models\SplineDyn\genMatrix\run\matrix_3p0_reduce1p0.nc','lon_matrix_beta');

plot(dist1,abs(beta1(1,:)))
hold on
plot(dist2,abs(beta2(1,:)))
hold on
plot(dist3,abs(beta3(1,:)))
hold on

eps = 0.8;

rbf_base = @(r,eps) exp(-(eps.*r).^2);
% rbf_base = @(r,eps) r.^eps;
% rbf_base = @(r,eps) sqrt(1+(eps.*r).^2);

rbf_base_deriv = @(r,eps,xd) eps^2.*exp(-(eps.*r).^2).*xd;
% rbf_base_deriv = @(r,eps) 2*eps*r.^(2*eps-1);
% rbf_base_deriv = @(r,eps) eps^2*r ./ sqrt(1+(eps*r).^2);

new_dist1 = bsxfun(@circshift,dist1,0:size(dist1,1)-1);
new_dist2 = bsxfun(@circshift,dist2,0:size(dist2,1)-1);
new_dist3 = bsxfun(@circshift,dist3,0:size(dist3,1)-1);

A1 = rbf_base(new_dist1,eps);
A2 = rbf_base(new_dist2,eps);
A3 = rbf_base(new_dist3,eps);

xd1 = dist1;
xd2 = dist2;
xd3 = dist3;

xd1(end/2+1:end) = -xd1(end/2+1:end);
xd2(end/2+1:end) = -xd2(end/2+1:end);
xd3(end/2+1:end) = -xd3(end/2+1:end);

B1 = rbf_base_deriv(dist1,eps,xd1);
B2 = rbf_base_deriv(dist2,eps,xd2);
B3 = rbf_base_deriv(dist3,eps,xd3);

D1 = A1\B1;
D2 = A2\B2;
D3 = A3\B3;

plot(dist1,abs(D1),'LineWidth',3)
hold on
plot(dist2,abs(D2),'LineWidth',3)
hold on
plot(dist3,abs(D3),'LineWidth',3)
hold on

f1 = sin(lon1);
f2 = sin(lon2);
f3 = sin(lon3);

fd1=D1.*f1;
fd2=D2.*f2;
fd3=D3.*f3;