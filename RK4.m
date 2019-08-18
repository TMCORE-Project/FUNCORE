% 4th order Runge-Kutta
function [stat_new] = RK4(stat_old,mesh,dt)

tend(1) = spatial_operator(stat_old,mesh);
stat(2) = update_stat(stat_old,tend(1),0.5*dt);

tend(2) = spatial_operator(stat(2),mesh);
stat(3) = update_stat(stat_old,tend(2),0.5*dt);

tend(3) = spatial_operator(stat(3),mesh);
stat(4) = update_stat(stat_old,tend(3),dt);

tend(4) = spatial_operator(stat(4),mesh);

tend_new.u  = ( tend(1).u  + 2*tend(2).u  + 2*tend(3).u  + tend(4).u ) / 6;
tend_new.v  = ( tend(1).v  + 2*tend(2).v  + 2*tend(3).v  + tend(4).v ) / 6;
tend_new.w  = ( tend(1).w  + 2*tend(2).w  + 2*tend(3).w  + tend(4).w ) / 6;
tend_new.gh = ( tend(1).gh + 2*tend(2).gh + 2*tend(3).gh + tend(4).gh) / 6;

stat_new = update_stat(stat_old,tend_new,dt);