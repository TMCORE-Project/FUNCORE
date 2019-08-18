function stat_new = update_stat(stat_old,tend,dt)
stat_new.u  = stat_old.u  + dt * tend.u;
stat_new.v  = stat_old.v  + dt * tend.v;
stat_new.w  = stat_old.w  + dt * tend.w;
stat_new.gh = stat_old.gh + dt * tend.gh;
