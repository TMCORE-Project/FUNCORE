% Evaluates the RHS (spatial derivatives) for the Cartesian RBF formulation of 
% the shallow water equations with projected gradients.
% This function applies to Test Case 5, which contains bottom topography
function tend = spatial_operator(stat,mesh)

H = [stat.u, stat.v, stat.w, stat.gh];

% Extract out some constants from the atm structure to make the code
% easier to read.
x = mesh.xCell;
y = mesh.yCell;
z = mesh.zCell;
f = mesh.f;
a = mesh.a;

% Compute the (projected) Cartesian derivatives applied to the velocity
% and geopotential.
parfor i = 1:4
Tx(:,i) = deriv(H(:,i),mesh.dx,mesh.kdtree,a);
Ty(:,i) = deriv(H(:,i),mesh.dy,mesh.kdtree,a);
Tz(:,i) = deriv(H(:,i),mesh.dz,mesh.kdtree,a);
end

%
% This is the computation for the right hand side of the (Cartesian) 
% momentum equations.
%
p = -(H(:,1).*Tx(:,1) + H(:,2).*Ty(:,1) + H(:,3).*Tz(:,1) + (f).*(y.*H(:,3) - z.*H(:,2)) + Tx(:,4));
q = -(H(:,1).*Tx(:,2) + H(:,2).*Ty(:,2) + H(:,3).*Tz(:,2) + (f).*(z.*H(:,1) - x.*H(:,3)) + Ty(:,4));
s = -(H(:,1).*Tx(:,3) + H(:,2).*Ty(:,3) + H(:,3).*Tz(:,3) + (f).*(x.*H(:,2) - y.*H(:,1)) + Tz(:,4));

% Project the momentum equations onto the surface of the sphere.
td(:,1) = mesh.p_u(:,1).*p + mesh.p_u(:,2).*q + mesh.p_u(:,3).*s;
td(:,2) = mesh.p_v(:,1).*p + mesh.p_v(:,2).*q + mesh.p_v(:,3).*s;
td(:,3) = mesh.p_w(:,1).*p + mesh.p_w(:,2).*q + mesh.p_w(:,3).*s;

% Right-hand side for the geopotential (Does not need to be projected, this
% has already been accounted for in the DPx, DPy, and DPz operators for
% this equation).
td(:,4) = -( H(:,1).*(Tx(:,4) - mesh.dghs(:,1)) ...
           + H(:,2).*(Ty(:,4) - mesh.dghs(:,2)) ...
           + H(:,3).*(Tz(:,4) - mesh.dghs(:,3)) ...
           + (H(:,4)-mesh.ghs).*(Tx(:,1) + Ty(:,2) + Tz(:,3)));

% Apply the hyper-viscosity, either once or twice.
parfor i = 1:4
    HH = H(:,i);
    td(:,i) = td(:,i) + sum(mesh.L .* HH(mesh.kdtree),2);
end
% F = F - L*(L*H);
% F = reshape(F,4*nCells,1);

tend.u  = td(:,1);
tend.v  = td(:,2);
tend.w  = td(:,3);
tend.gh = td(:,4);
