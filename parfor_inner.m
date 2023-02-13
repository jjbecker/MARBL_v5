function [sim, bgc, time_series, forcing, ierr, x , fnrm] = parfor_inner(sim, MTM, tracer_str)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     tracer_str = sim.tracer_loop(par_idx)
timer_total = tic;
tName = tracer_names(0);    % no CISO tracers

% ALWAYS punt "ALT" methods that do NOT depend on any other tracers...
% remove duplicates and make sure all choices are valid...
% sim.selection = [ find( strcmp(tName,tracer_str) ) ];
% sim.selection(ismember(sim.selection, [9,11]))=[];
% sim.selection = unique(sort(sim.selection));
sim.selection = unique(sort(find( strcmp(tName,tracer_str) )));
fprintf('%s.m: Selected tracer(s): #%d, "%s"\n', mfilename, sim.selection, string(tName(sim.selection)'));

% FIXME: endless headache of MARBL threads! Kill any leftover MEX running.
clear functions
clear calc_G calculate_forcing phi time_step_ann % clear "persistent" vars



fprintf('%s.m: Reading (Matlab) OFFLINE sim restart file with tracers and transports: %s\n', mfilename, sim.inputRestartFile);
if ~isfile(sim.inputRestartFile)
    error("missing file or typo in name of inputRestartFile")
end

% remember brsola() "sol" is x0 value. x1 value is NOT last col of
% x_hist; it is sol !!!
% First relax iteration of x0 gives same x1 as sol run.
% To be useful num_single_tracer_relax_iters >= 2 if using x0_sol, but OK for x1_sol
%%%%%% End of "inputs"

% Stuff below is not a simulation input. It is code to setup grids, etc a
% forward integration, or NK(), or...
%

%%%%%%
sim.runInParallel = 0;      % parfor can't use spmd inside, at least I can not make that work
[sim, bgc, ~, time_series, forcing] = init_sim(sim);
sim = calc_global_moles_and_means(bgc, sim);


%%%%%%

% =============== This is NK solver code ================

toc(timer_total)

fprintf('%s.m: Start Newton (Broyden Method) solver for tracer: ( %s ) @ %s\n',mfilename, string(tracer_str),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
timer_PQ_init_solve_relax_fwd = tic;
% fprintf('%s.m: Solving for tracer: %s\n', mfilename, string(tracer_str));


%%%%
% % DEBUG: add a known spike to see if it gets removed by Nsoli()spike_size = 1234.567;
% spike_size = 0;
% % spike_size = 1234.567;
% fprintf('\n**** %s.m: Adding a spike of %g to water parcel (%d,%d,%d) ****\n\n', mfilename, spike_size, sim.time_series_loc, sim.time_series_lvl,sim.selection);
% bgc.tracer(sim.time_series_loc,sim.time_series_lvl,sim.selection) = spike_size +bgc.tracer(sim.time_series_loc,sim.time_series_lvl,sim.selection);
%%%%

% confusing? YES
% tracer = [7881, 60, 32]
% x, sol = [379913]
% c      = [379913, 32]
% c0     = 12157216]

sz = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
c0 = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
c  = reshape(c0, sz);

% Solve only on selected tracers

x0 = c(:,sim.selection);    % initial condition for Nsoli()
x0 = x0(:);                 % unitless


if( sim.forwardIntegrationOnly )

    fprintf('%s.m: forward integration ONLY\n',mfilename);
    if sim.debug_disable_phi
        fprintf('\n\n\t%s.m: ********* phi() is short circuited: set x = 0*x0 *********\n\n',mfilename)
        x = 0 *c(:,1);
        fnrm = 0;
    end

else % Solve for selected tracer

    % NK always using 1 year integration

    if sim.phi_years ~= 1
        disp('ERROR: NK requires sim.phi_years = 1')
        keyboard
    end

    PQ_inv = calc_PQ_inv(sim, bgc, time_series, forcing, MTM);

    f = @(x) calc_G(x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv);
    [ierr, fnrm, myRestartFile_x0, x0_sol, c0, sim, bgc] = marbl_solve(x0, c0, sim, bgc, f);
    %     f0=feval(f,x0);
    %     [ierr, fnrm, myRestartFile_x0, x0_sol, c0, sim, bgc] = marbl_solve(x0, c0, sim, bgc, f, f0);


    x = x0_sol;     % FIXME: use x0 or x1 of marbl_solve?

    if sim.num_single_tracer_relax_iters > 0
        % keyboard
        [x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv, myRestartFile_relaxed] = ...
            marbl_relax(x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv);
    end % relax step

end % Solve for selected tracer

% be sure to shutdown MEX
mex_marbl_driver('shutdown');


elapsedTime_all_loc = toc(timer_PQ_init_solve_relax_fwd);
disp(' ');
% disp([mfilename,'.m: finished ', strjoin(tName(sim.selection))])
fprintf('%s.m: Finished Newton (Broyden Method) solver for tracer: ( %s ) @ %s\n',mfilename, strjoin(tName(sim.selection)),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
disp([mfilename,'.m: Runtime ', num2str(elapsedTime_all_loc, '%1.0f'),' (s) or ', num2str(elapsedTime_all_loc/60, '%1.1f'), ' (m)'])
%     disp(['Runtime per location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps/sim.domain.num_wet_loc*1000, '%1.2f'), ' (ms) MARBL, advection, diffusion, mfactor()'])
%     disp(['Runtime all location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps, '%1.2f'),                    ' (s)  MARBL, advection, diffusion, mfactor()'])
%     disp(['Runtime all location per sim year : ', num2str(elapsedTime_all_loc/60/1440/sim.tot_t*sim.const.sec_y, '%1.2f'), ' (d/y_sim)'])
disp(['Simulation speed: ', num2str(sim.tot_t/elapsedTime_all_loc/sim.const.days_y, '%1.1f'), ' (sim y/d) aka (SYPD)'])
disp(' ');

%%
if sim.debug_disable_phi
    fprintf('\n\n\t%s.m: ********* phi() is short circuited skip read of inputRestartFile  *********\n',mfilename);
    ierr = randi([1000 2000]);
    fprintf('\t%s.m: ********* phi() is short circuited force ierr to %d             *********\n\n',mfilename, ierr);
else
    sim.inputRestartFile = myRestartFile_x0;
end % if
end % end function