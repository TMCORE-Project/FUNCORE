function [stat_new] = temporal_integration(stat_old,mesh,dt,scheme)

if strcmp(scheme,'RK4')
    [stat_new] = RK4(stat_old,mesh,dt);
end

