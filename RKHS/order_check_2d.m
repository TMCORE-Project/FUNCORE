clc
clear

eps       = 5;
base_type = 1;
poly_deg  = 5;

res(1) = 0.25;
res(2) = 0.5;
res(3) = 1;
res(4) = 2;

d2r = pi/180;
r2d = 180/pi;

% Select base function
if base_type == 1
    rbf_base = @(r,eps) exp(-(eps.*r).^2);
    rbf_base_deriv = @(r,eps,xd) -2*eps^2.*exp(-(eps.*r).^2).*xd;
elseif base_type == 2
    rbf_base = @(r,eps) r.^eps;
    rbf_base_deriv = @(r,eps,xd) eps * r.^(eps-2).*xd;
elseif base_type == 3
    rbf_base = @(r,eps) sqrt(1+(eps.*r).^2);
    rbf_base_deriv = @(r,eps,xd) eps.^2 ./ sqrt(1+(eps*r).^2).*xd;
end

for ires = 1:size(res,2)
    x1d = -10:res(ires):10;
    y1d = -10:res(ires):10;
    
    nx = size(x1d,2);
    ny = size(y1d,2);
    
    [x2d,y2d] = meshgrid(x1d,y1d);
    
    x = reshape(x2d,[],1);
    y = reshape(y2d,[],1);
    
    x(1) = 0;
    y(1) = 0;
    
    w1 = gen_weights_2d         (x,y,eps,poly_deg,base_type);
    w2 = RBF_FD_PHS_pol_weights (x,y,eps,poly_deg);
    
    f   = 1 + sin(4*x) + cos(3*x) + sin(2*y);
    f2d = reshape(f,nx,ny);
    
    dfdx  = 4*cos(4*x) - 3*sin(3*x);
    dfdy  = 2*cos(2*y);
    
    dfdx1(ires) = sum(w1(:,1).*f);
    dfdx2(ires) = sum(w2(:,1).*f);
end

L2_old = abs(dfdx1-dfdx(1))/abs(dfdx(1));
L2_new = abs(dfdx2-dfdx(1))/abs(dfdx(1));
order_old  = log(L2_old(2)/L2_old(1))/log(2);
order_new  = log(L2_new(2)/L2_new(1))/log(2);