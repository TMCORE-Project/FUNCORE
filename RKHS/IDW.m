%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INVERSE DISTANCE WEIGHT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vint]=IDW(xc,yc,vc,x,y,eps,r1,r2,bf_opt)
%%%%%%%%%%%%%%%%%
%%% INPUTS
%xc = stations x coordinates (columns) [vector]
%yc = stations y coordinates (rows) [vector]
%vc = variable values on the point [xc yc]
%x = interpolation points  x coordinates [vector]
%y = interpolation points y coordinates [vector]
%e = distance weight
%r1 --- 'fr' = fixed radius ;  'ng' = neighbours
%r2 --- radius lenght if r1 == 'fr' / number of neighbours if  r1 =='ng'
%bf_opt --- base function option, 1 for r^eps, 2 for gaussian function exp(-(eps*r).^2)
%%% OUTPUTS
%Vint --- Matrix [length(y),length(x)] with interpolated  variable values
%%% EXAMPLES
%%% --> V_spa=IDW(x1,y1,v1,x,y,-2,'ng',length(x1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simone Fatichi -- simonef@dicea.unifi.it
%   Copyright 2009
%   $Date: 2009/06/19 $
%   $Updated 2012/02/24 $
%   Modified by Lilong Zhou Aug,19,2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vint=zeros(length(x),length(y));
% xc=reshape(xc,1,length(xc));
% yc=reshape(yc,1,length(yc));
% vc=reshape(vc,1,length(vc));
xs=reshape(xc,[],1);
ys=reshape(yc,[],1);
vs=reshape(vc,[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rbf
rbf_base = @(eps_rbf,r) exp(-(eps_rbf.*r).^2);
eps_rbf  = 1;
r_1d     = pdist([xs,ys]);
r        = squareform(r_1d);
phi      = rbf_base(eps_rbf,r);
c        = phi \ vs;

% vs_ext            = vs;
% vs_ext(end+1)     = 0;
% phi_ext           = phi;
% phi_ext(end+1,:)  = 1;
% phi_ext(:,end+1)  = 1;
% phi_ext(end,end)  = 0;
% c_ext             = phi_ext \ vs_ext;
% c                 = c_ext(1:end-1);

% rbf_val = rbf(xs, ys, vs, x, y, 'gaussian',0,-4);

if  strcmp(r1,'fr')
    if  (r2<=0)
        disp('Error: Radius must be positive')
        return
    end
    for j=1:length(y)
        for i=1:length(x)
            D = sqrt((x(i)-xs).^2 +(y(j)-ys).^2);
            vcc = vs(D<r2);
            D   = D(D<r2);
            cf  = c(D<r2);
            wV  = vcc.*bf(D,eps,bf_opt);
            w   = bf(D,eps,bf_opt);
            
            inf_pos = isinf(w);
            if sum(inf_pos)~=0
                wV = vcc(inf_pos);
                w  = 1;
            end
            
            if isempty(D)
                wV=NaN;
            else
                if bf_opt == 4
                    wV=sum(wV.*sum(cf.*rbf_base(eps_rbf,D)))/sum(w);
                else
                    wV=sum(wV)/sum(w);
                end
            end
            
            Vint(i,j)=wV;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    if (r2 > length(vs)) || (r2<1)
        disp('Error: Number of neighbours not congruent with data')
        return
    end
    for j=1:length(y)
        for i=1:length(x)
            D= sqrt((x(i)-xs).^2 +(y(j)-ys).^2);
            [D,I]=sort(D);
            vcc = vs(I);
            cf  = c(1:r2);
            wV  = vcc(1:r2).*bf(D(1:r2),eps,bf_opt);
            w   = bf(D(1:r2),eps,bf_opt);
            
            inf_pos = isinf(w);
            if sum(inf_pos)~=0
                wV = vcc(inf_pos);
                w  = 1;
            end
            
            if bf_opt == 4
                wV=sum(wV)/sum(w);
            else
                wV=sum(wV.*sum(cf.*rbf_base(eps_rbf,D)))/sum(w);
            end
            
            Vint(i,j)=wV;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
end

% base function
function f = bf(r,eps,bf_opt)
if bf_opt == 1 || bf_opt == 4
    if eps > 0
        disp('eps must be negative for inverse distance scheme')
        reture
    end
    
    f = r.^eps;
    
elseif bf_opt == 2
    if eps > 0
        disp('eps must be positive for gaussian scheme')
        reture
    end
    
    f= exp(-(eps*r).^2);
    
elseif bf_opt == 3
    if eps > 0
        disp('eps must be negative for inverse distance scheme')
        reture
    end
    
    dr = eps - r;
    dr(dr<0) = 0;
    
    f = ( dr ./ (eps*r) ).^2;
end
end