function [x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv, myRestartFile_relaxed] = marbl_relax(num_relax_iterations, x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fprintf('\n%s.m: Now we "relax" or forward integtration single variable solution a few iterations: %s\n',mfilename,datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
fprintf('%s.m: %d time steps, each time step is %1.1f (h), simulating %1.1f years\n', ...
    mfilename, sim.T/sim.dt, sim.dt/sim.const.sec_h, sim.tot_t /sim.const.sec_y)

% FIXME: keep history of relax steps???
x_histx = x;
r_hist = zeros(size(x));
it_histx = zeros(num_relax_iterations,1);

% % remember! First relax iteration of x0_sol" is same as x1 of sol, so in
% % that case, to be useful, num_relax_iterations >= 2.
% But if using x1_sol num_relax_iterations >= 1 is ok.

for itc = 1:num_relax_iterations

    fprintf("\n%s.m: starting relaxation year #%d of %d\n", mfilename, itc, num_relax_iterations)

    % Note x = x1, not "x0" which is normal thing for sol iterations
    % Note c0 is correct
    %                 [r, G, x1] = calc_G(x,c0,sim,bgc_sol_x1,time_series,forcing,MTM,PQ_inv);
    [r, G, x1] = calc_G(x,c0,sim,bgc,time_series,forcing,MTM,PQ_inv);

    % DEBUG
    fprintf('%s.m: itc %d norm(G) = %g\n',mfilename, itc, norm(G))
    fprintf('%s.m: itc %d norm(r) = %g\n',mfilename, itc, norm(r))
    x_histx = [x_histx,x];
    r_hist  = [r_hist, r];
    it_histx(itc,1) = norm(r);

    x = x1;

end % relax loop

% if we did NOT relax, then x = x1_sol, or x0_sol FIXME
%   If we relaxed x1, x = x2_sol, etc, etc


% Note c0 is correct
nsoli_relax = replaceSelectedTracers(sim, c0, x, sim.selection);
bgc_relax = bgc;
bgc_relax.tracer = nsoli2bgc(sim, bgc_relax, nsoli_relax);

bgc = bgc_relax;    % final answer or input to fwd integration...

tName = tracer_names(0);    % no CISO tracers
myRestartFile_relaxed = sprintf('%s/restart_%d_%s_relax_x1.mat', sim.outputRestartDir, round(sim.start_yr), strjoin(tName(sim.selection)));
[sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile_relaxed);

end