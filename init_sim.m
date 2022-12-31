function [sim, bgc, bgc_struct, time_series, forcing] = init_sim(sim)
%UNTITLED7 Summary of this function goes here

fprintf('Starting %s.m...\n',mfilename)
timer_total = tic;
tic;

% set size of OCEAN -not- world grids...  (num_wet_loc)
% Copy initial state, tracer, and so on -NOT- M3d

% Define time_step in seconds and number of time steps
% Note: for really long sim, like 1,000 years have to think carefully about
% sideral days and all that stuff. Here we are talking "mean solar day"
%
% https://www.timeanddate.com/time/earth-rotation.html
%
% "If every day were as long as today, a negative leap second would have
% to be added every 1398.99 days....

sim.const.sec_m = 60;
sim.const.sec_h = 60    *sim.const.sec_m;
sim.const.sec_d = 24    *sim.const.sec_h;   % mean solar "appearant" day
% Use Laskar's expression mean tropical year 2000
% sim.const.days_y = 365.242189; % mean tropical year = 4 seasons
% avoid leap year insanity for cyclostationary?
sim.const.days_y = 365;

sim.const.sec_y = sim.const.days_y *sim.const.sec_d;

sim.dt = sim.const.sec_h *sim.time_step_hr;
sim.num_time_steps = round ( sim.phi_years *sim.const.sec_y /sim.dt );
sim.T  = sim.num_time_steps * sim.dt;
sim.tot_t = sim.dt*sim.num_time_steps;  % used only for debug output


% % % Define gird dimensions in meters
[sim, bgc_struct] = calculate_depth_map_and_volumes(sim);
%%%%%%%%%%%%%%

% Initialize MARBL with "domain" aka grid dimensions.
% FIXME: initialize MARBL and read number of trancers, etc ??? forces
% dimensions here match, however if dimensions here do not match, then
% most of "names, units, etc" files in "utils" directory are going
% to be wrong too. Use dimension in utils files and then check that
% they match when MARBL is initialized.

load(sim.inputRestartFile,'tracer','state','forcing');
bgc.tracer       = tracer;      clear tracer;
% load(sim.inputRestartFile,'tracer','state','forcing','old_tracer');
% bgc.old_tracer   = old_tracer;  clear old_tracer;

if (sim.verbose_debug) 
    disp('Initializing global grids for tracer, tendency, etc...')
end

if (sim.runInParallel)
    sim.number_of_threads = 4; % only 4 on laptop or 12 on GP supported
else
    fprintf('\n\n%s.m: NOT running in parallel\n\n\n', mfilename);
end

disp('Starting client MEX...')
[sim, bgc_struct] = init_marbl(sim.marbl_file,sim, bgc_struct, forcing(1).surf_forcing);
disp('...client MEX is running')

% Make sure size of tracers, e.g. CISO, matches sim we are about
% to run, not size we used to make default (which is always
% with CISO)
%
% This is about creating correct sizes of array for whole
% ocean; its not about loading in values from forcing file...

% We are given actual values of bgc.tracer
% bgc.tracer    = makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.tracer);

bgc.kmt         = sim.domain.bottom_lvl(sim.domain.wet_loc);
bgc.forcing     = makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.forcing);
bgc.surf_forcing= makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.surf_forcing);
% FIXME: surface_state???
if exist('state','var')
    %     bgc.surf_state  = surf_state;
    bgc.state       = state;
else
    bgc.state       = makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.state);
end
bgc.tendency    = makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.tendency);
% bgc.sfo         = makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.sfo);
% bgc.diag        = makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.diag);
% bgc.surf_diag   = makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.surf_diag);
% bgc.surf_flux   = makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.surf_flux);
% river flux is not given to MARBL, only used in time_step()
bgc.river_flux  = makeLinearGrid(sim.domain.wet_loc, bgc_struct.size.grd, bgc_struct.value.river_flux);

% bgc.accumulate  = 0 *packMarbl( bgc.tracer, sim.domain.iwet_JJ );


% Limit default value? If they are that far off, why is it a default?
% FIXME:% Interpolation can produce small negative oscillations.
% bgc.tracer(:)   = max(1e-15, bgc.tracer(:));

sim.bgc_struct_base = bgc_struct;
time_series = init_time_series(sim, bgc, bgc_struct);

toc
if (sim.verbose_debug) 
    % FIXME: super useful check of forcing flux for Fe and so on, but super
    % verbose for sure.
    % checkRestartFile(sim, bgc, forcing)
end

%%
if (sim.runInParallel)
    tic;

    % First shut down client MARBL previously used to get dim of MEX arrays
    disp('Shutting down client MEX (aka serial) because it interfers with parallel MEX threads...')
    mex_marbl_driver('shutdown');
    % LOCAL_MARBL_IS_RUNNING = 0;
    disp('...client MEX shutdown inside if (sim.runInParallel) worked.')

    % Also shutdown existing worker pool, if any

    disp('Shutting down existing pool of threads, if any...')
    delete(gcp('nocreate'));
    sim.number_of_threads = min(sim.number_of_threads, sim.domain.num_wet_loc);    % >= num of rows or crash

    %     parpool('threads')
    parpool(sim.number_of_threads);
    ticBytes(gcp)
    toc
    %     interior_base_client = interior_base;
    if (sim.runInParallel)
        % FIXME: this is not place for this; not needed; etc
        % synchronize all labs
        %    labBarrier;
        mpiprofile reset
        % mpiprofile on
    end

    spmd(sim.number_of_threads)
        mex_marbl_driver('read_settings_file', sim.marbl_file);

        if sim.lciso_on
            mex_marbl_driver('put setting', 'ciso_on = .true.');
        end

        mex_marbl_driver ( 'init', ...
            sim.domain.dzt *sim.domain.MARBL_depth_per_m, ...
            sim.domain.zw  *sim.domain.MARBL_depth_per_m, ...
            sim.domain.zt  *sim.domain.MARBL_depth_per_m );
    end
    elapsedTime = toc;
    disp(['start threads, boot Matlab on them: ', num2str(elapsedTime, '%1.1f'),' (s)']);
end

fprintf('%s.m finished: Total time to read init files, start parallel threads, etc : %s sec\n', mfilename,num2str(toc(timer_total),'%.1f'))

end % init_sim.m
