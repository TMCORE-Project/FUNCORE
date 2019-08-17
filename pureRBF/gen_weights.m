function w = gen_weights (x,y,z,eps)
% Input parameters
% x,y,z Column vectors with stencil node locations; approximation to
% be accurate at x(1),y(1),z(1)
% m Power of r in RBF fi(r) = r^m, with m odd, >= 3.
% d Degree of supplementary polynomials (d = -1 no polynomials)
%
% Output parameter
% w Matrix with three columns, containing weights for d/dlon, d/dlat

xd = x-x(1);
yd = y-y(1);
zd = z-z(1);% Shift nodes so stencil centered at origin

% ------ RBF part --------------------------------------------------------
coord  = [x,y,z];
r_1d   = pdist(coord);
r      = squareform(r_1d);

phi    = exp( - eps^2 * r.^2 );

dphidx = -2 * eps^2 .* phi(:,1) .* xd;
dphidy = -2 * eps^2 .* phi(:,1) .* yd;
dphidz = -2 * eps^2 .* phi(:,1) .* zd;

L      = [dphidx,dphidy,dphidz]; % RHSs

w      = phi\L;
