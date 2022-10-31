function [sim, bgc, bgc_struct, time_series, forcing] = init_sim(marbl_file, restart_file, sim, phi_years, dt)
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

sim.dt = sim.const.sec_h *dt;
sim.num_time_steps = round ( phi_years *sim.const.sec_y /sim.dt );
sim.T  = sim.num_time_steps * sim.dt;

% Define gird dimensions in meters

sim.domain.zt   = sim.grd.zt;
sim.domain.zw   = sim.grd.zw+sim.grd.dzt; % in MARBL zw is bottom, not top...
sim.domain.dzt  = sim.grd.dzt;
sim.domain.dzw  = sim.grd.dzw;

% ---> MARBL uses units of - CM - for depths, Matlab code uses meters!
sim.domain.MARBL_depth_per_m = 100;     % 1 meter here is 100 cm MARBL uses

% Create a struct that has meta data like size, names, and units.

bgc_struct= init_bgc_struct(sim);


[sim.domain.iwet_JJ, sim.domain.wet_loc, sim.domain.iwet_FP, sim.domain.num_wet_loc, ...
    sim.domain.bottom_lvl, sim.domain.iColLvlFromFp] = createIwetEtc(sim.domain.M3d);

% 2d matrix of "depth", e.g. 0:10900 meter <-- needed by MARBL
% FIXME: is this actually used anywhere???
% Only thing tricky is that Matlab array cant have index of 0...

zw = cumsum(sim.grd.dzt')';     % bottom of grid, not depth at given local

mapping = [1:length(zw); zw]';
mapping = [0 0; mapping];           % add row for "dry" with indx == 0
mapping(:, 1) = mapping(:, 1) +1;   % +1 for row with "dry" value
LUT(mapping(:, 1)) = mapping(:, 2);
sim.domain.bottom_depth = LUT(sim.domain.bottom_lvl+1);   % +1 for row with "dry" value

sim.domain.dVt    = sim.grd.dV;
sim.domain.dVt    = sim.domain.M3d .* sim.domain.dVt;
sim.domain.V      =         sum(sim.domain.dVt, 'all') ;
% with depth: layers get thicker, but area decreases. Volume per level
% almost constant but sort of random
sim.domain.V_lvl  = squeeze(sum(sim.domain.dVt, [1 2])); % vol at given iLvl
sim.domain.V_lat  = squeeze(sum(sim.domain.dVt, [3 2])); % vol at given iLat
sim.domain.V_lon  = squeeze(sum(sim.domain.dVt, [3 1])); % vol at given iLon

sim.domain.dVt_FP = sim.domain.dVt(sim.domain.iwet_FP);

%%%%%%%%%%%%%%

% Initialize MARBL with "domain" aka grid dimensions.
% FIXME: initialize MARBL and read number of trancers, etc ??? forces
% dimensions here match, however if dimensions here do not match, then
% most of "names, units, etc" files in "utils" directory are going
% to be wrong too. Use dimension in utils files and then check that
% they match when MARBL is initialized.

load(restart_file,'tracer','state','forcing');
bgc.tracer       = tracer;      clear tracer;
% load(restart_file,'tracer','state','forcing','old_tracer');
% bgc.old_tracer   = old_tracer;  clear old_tracer;

if (sim.verbose_debug) disp('Initializing global grids for tracer, tendency, etc...'), end;
[sim, bgc_struct] = init_marbl(marbl_file,sim, bgc_struct, forcing(1).surf_forcing);

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
time_series = init_time_series(sim, bgc_struct);

toc
if (sim.verbose_debug) checkRestartFile(sim, bgc, forcing), end
%%
if (sim.runInParallel)
    tic;

    % First shut down client MARBL previously used to get dim of MEX arrays
    disp('Shutting down client MEX (aka serial) because it interfers with parallel MEX threads')
    mex_marbl_driver('shutdown');

    % Also shutdown existing worker pool, if any

    disp('Shutting down existing pool, if any...')
    delete(gcp('nocreate'));
    sim.number_of_threads = min(sim.number_of_threads, sim.domain.num_wet_loc);    % >= num of rows or crash

    %     parpool('threads')
    parpool(sim.number_of_threads);
    ticBytes(gcp)
    toc
    %     interior_base_client = interior_base;
    if (sim.runInParallel)
        % sim.number_of_threads = 4;
        % FIXME: this is not place for this; not needed; etc
        % synchronize all labs
        %    labBarrier;
        mpiprofile reset
        % mpiprofile on
    end

    spmd(sim.number_of_threads)
        mex_marbl_driver('read_settings_file', marbl_file);

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
