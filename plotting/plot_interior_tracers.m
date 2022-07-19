function myFig = plot_interior_tracers(myFig, my_time, plot_layer, idx, sim, time_series, n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

dx = sim.const.sec_d/sim.dt;

t = (0:n-1) /dx;

plot_depth = sim.domain.zt(plot_layer);

% my_step = round(max( 1, min(round(my_time/sim.dt),sim.dt)))
% my_step = round(max( 1, min(round(my_time/sim.dt),sim.num_time_steps)))
my_step =       max( 1,     round(my_time/sim.dt));

nameUnits = strcat(sim.bgc_struct_base.name.tracer," ", sim.bgc_struct_base.unit.tracer);

% plot data v. time

myTitle = sprintf('Interior Tracers v. Time(d) @row %d, level %d, depth = %d(m)', ...
    sim.time_series_loc, plot_layer, round(plot_depth));
myData = squeeze(time_series.tracer(plot_layer, idx, 1:n))';
myFig = plot_log(myFig, myTitle, t, myData, nameUnits(idx), idx, false);

% plot data v. depth

myTitle = sprintf('Interior Tracers v. Depth (m) @row %d, iteration = %d, day %G', ...
    sim.time_series_loc, my_step, round(my_time/sim.const.sec_d,2));
myData = squeeze(time_series.tracer(:, idx, my_step));
myDepth = sim.domain.zt;
myFig = plot_log(myFig, myTitle, myDepth, myData, nameUnits(idx), idx, true);

end
