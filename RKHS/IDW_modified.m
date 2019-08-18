%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INVERSE DISTANCE WEIGHT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vint]=IDW(xc,yc,vc,x,y,e,r1,r2)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vint=zeros(length(x),length(y));
% xc=reshape(xc,1,length(xc));
% yc=reshape(yc,1,length(yc));
% vc=reshape(vc,1,length(vc));
xs=reshape(xc,[],1);
ys=reshape(yc,[],1);
vs=reshape(vc,[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            wV  = vcc.*(D.^e);
            w   = D.^e;
            if min(D)==0
                wV = vcc(D==0);
                w  = 1;
            end
            if isempty(D)
                wV=NaN;
            else
                wV=sum(wV)/sum(w);
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
            if min(D)==0
                wV = vcc(D==0);
                w  = 1;
            end
            [D,I]=sort(D);
            vcc = vs(I);
            wV  = vcc(1:r2).*(D(1:r2).^e);
            w   = D(1:r2).^e;
            wV  = sum(wV)/sum(w);
            Vint(i,j)=wV;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return