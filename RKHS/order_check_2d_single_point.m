clc
clear

eps       = 8;
base_type = 2;
poly_deg  = 45;
nSamples  = 49;

minRes    = 1;
nRes      = 2;

res(1) = minRes;
for ires = 2:nRes
    res(ires) = 2*res(ires-1);
end

d2r = pi/180;
r2d = 180/pi;

for ires = 1:nRes
    x1d = -10:res(ires):10;
    y1d = -10:res(ires):10;
    
    nx = size(x1d,2);
    ny = size(y1d,2);
    np = nx*ny;
    
    [x2d,y2d] = meshgrid(x1d,y1d);
    
    x = reshape(x2d,[],1);
    y = reshape(y2d,[],1);
    
    pts_i = [x,y];
    
    f   = 1 + sin(4*x) + cos(3*x) + sin(2*y);
    f2d = reshape(f,nx,ny);
    
    dfdx  = 4*cos(4*x) - 3*sin(3*x);
    dfdy  = 2*cos(2*y);
    
    dfdx_2d = reshape(dfdx,nx,ny);
    dfdy_2d = reshape(dfdy,nx,ny);
    
%     dfdx_n = zeros(np,1);
%     dfdy_n = zeros(np,1);
    %     parfor i = 1:np
    i = find(x==0&y==0);
    disp(['Calculate derivative weights for point # ',num2str(i)])
    pts_c = [x(i),y(i)];
    xc    = x(i);
    yc    = y(i);
    
    dist     = sqrt((x-xc).^2 + (y-yc).^2);
    
    [w,idx]  = gen_weights_2d (pts_i,pts_c,eps,poly_deg,nSamples,base_type);
    
    dfdx_n = sum(w(:,1).*f(idx));
    dfdy_n = sum(w(:,2).*f(idx));
    
    %     end
%     dfdx_n_2d = reshape(dfdx_n,nx,ny);
%     dfdy_n_2d = reshape(dfdy_n,nx,ny);
    
    L2_x(ires) = abs(dfdx_n - dfdx(i)) ./ abs(dfdx(i).^2);
    L2_y(ires) = abs(dfdx_n - dfdx(i)) ./ abs(dfdx(i).^2);
    
    if ires >= 2
        order_x(ires-1) = log(L2_x(ires)/L2_x(ires-1))/log(2);
        order_y(ires-1) = log(L2_y(ires)/L2_y(ires-1))/log(2);
    end
end
