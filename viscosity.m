function [ind_i, ind_j, weightsDx, weightsDy, weightsDz, weightsL] = viscosity(coord,tree,eps,fdsize,order,coef,dim,a)
%%% [D,L] = rbfmatrix_fd_tree(x,ep,alpha,fdsize,order,dim)
% IN:
% coord - nodes' coordinate
% eps - shape parameter
% fdsize - stencil size (good choices: 31, 50, 74, 101)
% order - L = L^order
% dim - dimension of Laplacian formula
% a - earth radius
% coef - diffusion coef for laplacian
% OUT:
% DPx - sparse differentiation matrix
% DPy - sparse differentiation matrix
% DPy - sparse differentiation matrix
% L - sparse dissipation matrix

N      = length(coord);

rbf = @(ep,rd2) exp(-ep^2*rd2);
drbf = @(ep,rd2) -2*ep^2*exp(-ep^2*rd2);

weightsDx = zeros(N,fdsize);
weightsDy = zeros(N,fdsize);
weightsDz = zeros(N,fdsize);
weightsL  = zeros(N,fdsize);

ind_i = zeros(N*fdsize,1);
ind_j = zeros(N*fdsize,1);

A = ones(fdsize,fdsize);
B = zeros(fdsize,4);

c = coef*eps^(2*order);

for j=1:N
    imat = tree(j,:)';
    outidx=(j-1)*fdsize+1:j*fdsize;
    
    dp = (coord(j,1)*coord(imat,1) + coord(j,2)*coord(imat,2) + coord(j,3)*coord(imat,3));
    
    rd2 = max(0,2*(1-coord(imat,1)*coord(imat,1).'-...
        coord(imat,2)*coord(imat,2).'-coord(imat,3)*coord(imat,3).'));
    rd2v = rd2(:,1);
    
    A(1:fdsize,1:fdsize) = rbf(eps,rd2);
    dr = drbf(eps,rd2v);
    B(1:fdsize,1:3) = [(coord(j,1)*dp-coord(imat,1)).*dr, ...
        (coord(j,2)*dp-coord(imat,2)).*dr, (coord(j,3)*dp-coord(imat,3)).*dr]/a;
    B(1:fdsize,4) = c*hyper(eps^2*rd2v,dim,order).*exp(-eps^2*rd2v);
    weights = A\B;

    t = sortrows([imat weights(1:fdsize,:)], 1);

    ind_i(outidx) = j;
    ind_j(outidx) = t(:,1);

    weightsDx(j,:) = t(:,2);
    weightsDy(j,:) = t(:,3);
    weightsDz(j,:) = t(:,4);
    weightsL (j,:) = t(:,5);
end

function p=hyper(ep2r2,d,k)
%%% laplacian to the power k of dimension d
% ep2r2 - ep^2*r^2
n = length(ep2r2);
P = zeros(n,k+1);
P(:,1) = 1;
P(:,2) = 4*ep2r2-2*d;
for j=3:k+1
    P(:,j) = 4*(ep2r2-2*j-d/2+4).*P(:,j-1) - ...
            8*(j-2)*(2*j+d-6)*P(:,j-2);
end
p = P(:,k+1);
