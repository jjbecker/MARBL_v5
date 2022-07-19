function myFig = plot_interior_diags(myFig, my_time, plot_layer, idx, name, sim,time_series, n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if max(idx)>size(time_series.diag,2)
    disp (' '); 
    idx_bad = idx;
    idx = idx(idx<size(time_series.diag,2));
    disp( ['non-CISO sim (???) called plot_interior_diags.m with idx that contained CISO diags: ', num2str(setdiff(idx_bad, idx), '%d ')]); 
    disp (' '); 
end

dx      = sim.const.sec_d/sim.dt;
t = (0:n-1) /dx;

plot_depth = sim.domain.zt(plot_layer);

my_step = max( 1, min(round(my_time/sim.dt),sim.num_time_steps));

% plot data v. time

myTitle = sprintf(append(name,' Interior Diags v. Time(d) @row %d, level %d, depth = %d(m)'), ...
    sim.time_series_loc, plot_layer, round(plot_depth));
myData = squeeze(time_series.diag(plot_layer, idx, 1:n))';
plot_log(myFig, myTitle, t, myData, sim.bgc_struct_base.name.diag(idx), idx, false);

myFig = myFig +1;

% plot data v. depth

myTitle = sprintf(append(name,' Interior Diags v. Depth (m) @row %d, iteration = %d, day %G'),...
    sim.time_series_loc, my_step, round(my_time/sim.const.sec_d,2));
myData = squeeze(time_series.diag(:, idx, my_step));
myDepth = sim.domain.zt;
plot_log(myFig, myTitle, myDepth, myData, sim.bgc_struct_base.name.diag(idx), idx, true);

myFig = myFig +1;

end
