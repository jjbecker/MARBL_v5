function [sim] = parfor_inner(sim, MTM, tracer_str)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     tracer_str = sim.tracer_loop(par_idx)
timer_total = tic;
tName = tracer_names(0);    % no CISO tracers

    % ALWAYS punt "ALT" methods that do NOT depend on any other tracers...
    % remove duplicates and make sure all choices are valid...
    sim.selection = [ find( strcmp(tName,tracer_str) ) ];
    sim.selection(ismember(sim.selection, [9,11]))=[];
    sim.selection = unique(sort(sim.selection));
    cstr = tName(sim.selection)';
    fprintf('%s.m: Selected tracer(s): #%d, "%s"\n', mfilename, sim.selection, string(cstr));

    % FIXME: endless headache of MARBL threads! Need to kill any
    % leftover MEX running on threads. This coincidentally clears persistent
    % variables in G() and phi().
    clear functions
    clear calc_G phi       % clear debugging counters usedin G() and phi()



    fprintf('%s.m: Reading (Matlab) OFFLINE sim restart file with tracers and transports: %s\n', mfilename, sim.inputRestartFile);
    if ~isfile(sim.inputRestartFile)
        error("missing file or typo in name of inputRestartFile")
    end

    % remember brsola() "sol" is x0 value. x1 value is NOT last col of
    % x_hist; it is sol !!!
    % First relax iteration of x0 gives same x1 as sol run.
    % To be useful num_relax_iterations >= 2 if using x0_sol, but OK for x1_sol
    %%%%%% End of "inputs"

    % Stuff below is not a simulation input. It is code to setup grids, etc a
    % forward integration, or NK(), or...
    %

    %%%%%%
sim.runInParallel = 0;      % parallel is hard to debug, but 2x faster
    [sim, bgc, ~, time_series, forcing] = init_sim(sim);

    % DEBUG: run on a few locations, rather global, MUCH faster debug something
    % sim.domain.num_wet_loc = 1; % comment out too run entire world

    % scale units to unity, using values of global mean of initial
    % need a column vector

    sim = calc_global_moles_and_means(bgc, sim);


    %%%%%%

    % =============== This is NK solver code ================

    toc(timer_total)

    fprintf('\n%s.m: Start Newton (Broyden Method) solver: %s\n',mfilename,datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
    timer_PQ_init_solve_relax_fwd = tic;
    fprintf('%s.m: Solving for tracer: %s\n', mfilename, string(tracer_str));


    % NK always using 1 year integration

    if sim.phi_years ~= 1
        disp('ERROR: NK requires sim.phi_years = 1')
        keyboard
    end

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

    moles_0 = global_moles(bgc.tracer, sim)';
    sz = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
    c0 = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
    c  = reshape(c0, sz);

    % Solve only on selected tracers

    x0 = c(:,sim.selection);    % initial condition for Nsoli()
    x0 = x0(:);                 % unitless


    if( sim.forwardIntegrationOnly )
        fprintf('%s.m: forward integration ONLY\n',mfilename);
    else
        % Solve for selected tracer
        if sim.recalculate_PQ_inv && ~(numel(sim.disabledPreconditoners)>0 && ismember(tName(sim.selection), sim.disabledPreconditoners))
            PQ_inv = calc_PQ_inv(sim, bgc, time_series, forcing, MTM);
        else
            tStart = tic;
            if sim.disable_ALL_Preconditioner || (numel(sim.disabledPreconditoners)>0 && ismember(tName(sim.selection), sim.disabledPreconditoners))
                fprintf('\n\n\t%s.m: ********* Replace preconditioner with 1 *********\n\n',mfilename)
                PQ_inv = 1
            else
                fprintf('\n%s.m: Loading ~30 GB(!) mfactored preconditioner PQ_inv from %s solution...\n', mfilename, strcat(string(tName(sim.selection))))
                PQ_inv = load (strcat(myDataDir(),'sol/',strjoin(tName(sim.selection)),'_QJ'), 'PQ_inv')
            end
            fprintf('%s.m: %1.0f (s) to init sim and load PQinv \n',mfilename, toc(tStart));
        end % calculate or load PQ_inv


        f = @(x) calc_G(x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv);
        f0=feval(f,x0);


        [ierr, myRestartFile_x0, x0_sol, c0, sim, bgc, time_series, forcing, MTM, PQ_inv] = ...
            marbl_solve(x0, c0, sim, bgc, time_series, forcing, MTM, PQ_inv, f, f0);



        % FIXME: use x0 or x1 of solution?
        x = x0_sol;

        if sim.num_relax_iterations > 0
            [x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv, myRestartFile_relaxed] = ...
                marbl_relax(x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv);
        end % relax step

    end % Solve for selected tracer

    % Next! allow ALL tracers, not just selection, to "relax' to solution.
    %    do pure forward integration for a while...

    for fwd_itc = 1:sim.num_forward_years

        fprintf("\n%s.m: starting forward integrate year #%d of %d\n", mfilename, fwd_itc, sim.num_forward_years)

        [sim, bgc, time_series] = phi(sim, bgc, time_series, forcing, MTM);

        sim.start_yr  = sim.start_yr+1;
        if mod(round(sim.start_yr), sim.yearsBetweenRestartFiles) == 0    % This runs after last time step of every 10 y

            myRestartFile_fwd = sprintf('%s/restart_%d_%s_fwd_x1.mat', sim.outputRestartDir, round(sim.start_yr),strjoin(tName(sim.selection)));
            [sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile_fwd);

        end
    end % fwd loop

    % this is my final answer!
    % always save my final answer
    myRestartFile_fwd = sprintf('%s/restart_%d_%s_fwd_x1.mat', sim.outputRestartDir, round(sim.start_yr),strjoin(tName(sim.selection)));
    [sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile_fwd);



    elapsedTime_all_loc = toc(timer_PQ_init_solve_relax_fwd);
    disp(' ');
    disp([mfilename,'.m: finished ', strjoin(tName(sim.selection))])
    disp(['Runtime: ', num2str(elapsedTime_all_loc, '%1.0f'),' (s) or ', num2str(elapsedTime_all_loc/60, '%1.1f'), ' (m)'])
    %     disp(['Runtime per location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps/sim.domain.num_wet_loc*1000, '%1.2f'), ' (ms) MARBL, advection, diffusion, mfactor()'])
    %     disp(['Runtime all location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps, '%1.2f'),                    ' (s)  MARBL, advection, diffusion, mfactor()'])
    %     disp(['Runtime all location per sim year : ', num2str(elapsedTime_all_loc/60/1440/sim.tot_t*sim.const.sec_y, '%1.2f'), ' (d/y_sim)'])
    disp(['Simulation speed: ', num2str(sim.tot_t/elapsedTime_all_loc/sim.const.days_y, '%1.1f'), ' (sim y/d) aka (SYPD)'])
    disp(' ');

    if sim.debug_disable_phi
        fprintf('\n\n\t%s.m: ********* phi() is short circuited skip inputRestartFile read  *********\n\n',mfilename)
    else
        sim.inputRestartFile = myRestartFile_x0;
    end
end