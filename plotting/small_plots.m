function small_plots (sim, time_series, n, row, lvl )

% As "learning tool" plot routine makes literally 100's plots.
% Edit "small_plots.m" code to make plots that are actually of interest.

tic;

nameUnits   = strcat(sim.bgc_struct_base.name.tracer," ", sim.bgc_struct_base.unit.tracer);
globalUnits = strcat(sim.bgc_struct_base.name.tracer," ", global_tracer_units(0));
surfaceDiagName  = sim.bgc_struct_base.name.surf_diag;


dt      = sim.dt;
tot_t   = dt*n;
my_time = tot_t;

% FIXME: do NOT pick midnight (aka 0 second of day) if you graph PAR
% my_time = tot_t -0.4 *sim.const.sec_d;

dx = sim.const.sec_d/dt;

t = (0:n) /dx;  % include initial tracer at 0 and updated tracers at n, so n+1 total points.

% surface tracers

fig = 1;
idx = 1:size(nameUnits');        % plot -ALL- tracers???
myTitle = sprintf('Surface Tracers v. Time (d) @row %d', row);
myData = squeeze(time_series.tracer(1,:,1:n+1))';
fig = plot_log(fig, myTitle, t, myData, nameUnits(idx), idx, false);

% surface diags

idx = 1:size(surfaceDiagName');
myTitle = sprintf('Surface Diags v. Time (d) @row %d', row);

if (sim.logDiags)
    myData = time_series.surf_diag(:,1:n+1)';        % plot -ALL- diags???
    fig = plot_log(fig, myTitle, t, myData, surfaceDiagName(idx), idx, false);
end

% Plot interior tracer time series

plot_layer = min(size(sim.domain.zt,2), lvl);

idx = 1:size(nameUnits');
fig = plot_interior_tracers(fig, tot_t, plot_layer, idx, sim, time_series,n+1); % include initial tracer at 0 and updated tracers at n, so n+1 total points.

% Plot surface observations

myTitle = 'SFO v. Time(d)';
idx = 1:4;
% TRICKY: we have last value of tracers but not SFO
t = (0:n-1) /dx;
myData = time_series.sfo(idx,1:n)';  
myName = sim.bgc_struct_base.name.sfo()';
fig = plot_log(fig, myTitle, t, myData, myName,idx, false);

% Plot interior accumualte

% idx = 1:size(nameUnits');
% fig = plot_accumulate(fig, tot_t, plot_layer, idx, sim, time_series,n);

% Plot interior diags

if (sim.logDiags)
    myTitle = 'Sample of';
%     idx = [ 5, 179, 149, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 151, 207, 208, 209, 210, 211 ];
    idx = [ 19, 31, 165, 166, 167, 168, 169, 170, 171, 172, 210, 211, 212, 213, 214];
    fig = plot_interior_diags(fig, my_time, plot_layer, idx, myTitle, sim, time_series,n);
end

% fig = 5;
idx = 1:size(nameUnits');        % plot -ALL- tracers???
myTitle = sprintf('***MOLES*** Global Volume Integrated Tracers(mole) v. Time (d)');
t = (0:n) /dx; % include initial tracer at 0 and updated tracers at n, so n+1 total points.
myData = squeeze(time_series.moles(:,1:n+1))'; % include initial tracer at 0 and updated tracers at n, so n+1 total points.
fig = plot_log(fig, myTitle, t, myData, globalUnits(idx), idx, false);

% fig = 6;
idx = 1:size(nameUnits');        % plot -ALL- tracers???
myTitle = sprintf('Global Volume Integrated(MARBL+River Tendency) ./Tracer *1y v. Time (d)');
% TRICKY: we have last value of tracers but not SFO
t = (0:n-1) /dx;
myData = sim.const.sec_y *(squeeze(time_series.Dmoles(:,1:n)) ./ squeeze(eps+time_series.moles(:,1:n)))';
fig = plot_log(fig, myTitle, t, myData, globalUnits(idx), idx, false);

% % fig = 9;
% idx = 1:size(nameUnits');        % plot -ALL- tracers???
% myTitle = sprintf('Fractional Global Volume Integrated (D(Tracer)/dt /(Tracer/dt) v. Time (d)');
% myData = squeeze(time_series.Dmoles(:,1:n))' ./ squeeze(time_series.moles(:,1:n))';
% fig = plot_log(fig, myTitle, t, myData, nameUnits(idx), idx, false);

return

%
%
%
%
% idx = [ 30, 23, 24, 338, 340, 179, 180, 174 ];
% fig = plot_interior_diags(fig, my_time, plot_layer, idx, 'POC', sim, time_series);
%
%
% idx = 90:116;
% fig = plot_interior_diags(fig, my_time, plot_layer, idx, 'Diaz Diags #1', sim, time_series);
%
% idx = 249:265;
% fig = plot_interior_diags(fig, my_time, plot_layer, idx, 'Diaz Diags #2', sim, time_series);
%
% idx = [ 327:334 369:372 385:386 ];
% fig = plot_interior_diags(fig, my_time, plot_layer, idx, 'Diaz Diags #3', sim, time_series);
%
% % plot_layer = 1;
% idx = [  277:288,301,303,305,306,309,311,313,316,319,321,324,327,329,332,337,338,373,374,377,379,381,383,385];
% fig = plot_interior_diags(fig, my_time, plot_layer, idx, '13C', sim, time_series);
%
% myTitle = 'Glitch in 13CO2 flux when pCO2 = atm (d)';
% idx = [ 13, 28:31, 37:40];
% myData = time_series.surf_diag(idx, :,);
% fig = plot_log(fig, myTitle, t, myData, surfaceDiagName(idx), idx, false);
%
% % Make a 3D plot of a small amount of time, say 7 days, starting day 10...
% % n_days = min( floor(sim.num_time_steps*dt/sim.const.sec_d), 7);
% n_days = min( floor(n*dt/sim.const.sec_d), 7);
% n_cnt = floor(n_days* sim.const.sec_d/dt);
% n_start = round(min(dt, 10* sim.const.sec_d)/dt);
% n_range = n_start:n_start +n_cnt -1;
%
% % -DIAG- #151 = PAR_avg = PAR at depth
% % great way to prove depths are in cm
% small_data = time_series.diag(:,:, n_range);
% idx = 151;
% fig = plot3dDiagTimeSeries  ( fig, small_data, sim, idx, n_cnt, dt);
%
% % Interior tracer #18 zoo plankton C
% % Great way to show Zoo do -NOT- move vertically during day/night
% small_data = time_series.tracer(:,:,n_range);
% idx = 18;
% fig = plot3dTracerTimeSeries( fig, small_data, sim, idx, n_cnt, dt);
%
% % Interior tracer #46 diaz C14 does show diurnal cycle
% idx = 46;
% fig = plot3dTracerTimeSeries( fig, small_data, sim, idx, n_cnt, dt);
%
% % Interior tracer #4
% idx = 4;
% fig = plot3dTracerTimeSeries( fig, small_data, sim, idx, n_cnt, dt);
% % "" but legible...
% figure(fig);
% foo = small_data(:,idx,plot_layer);
% t_unit = sim.const.sec_d;n_decimate = 1;
% t = decimate ( ( n_range ) *dt/t_unit, n_decimate );
% tracer_name  = sim.bgc_struct_base.name.tracer(idx);
% plot(t,foo);
% title("Tracer #"+idx+" "+tracer_name+" "+'plot_layer #'+plot_layer+ ...
%     ' depth '+sim.domain.zt(plot_layer)+" (m)",'Interpreter', 'none');
% xlabel('time (days)')
% ylabel(sim.bgc_struct_base.unit.tracer(idx))
% fig = 1+(fig);
%
%
% elapsedTime = toc;
% disp(['small_plots.m runtime: ', num2str(elapsedTime, '%1.0f'), ' (s)'])

end % small_plots.m
