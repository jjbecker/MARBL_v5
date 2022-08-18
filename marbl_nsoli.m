% "main" of cyclostationary transport version of MARBL,
%
%   matlab -nodisplay -nodesktop -nosplash -noFigureWindows -logfile batch.txt < marbl_nsoli.m &
%
% MARBL code being used is https://github.com/marbl-ecosys/MARBL
%
% function marbl()
% "main" of cyclostationary transport version of MARBL,

fprintf('Start of %s.m: %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));

% Need this to clear the "persistent" variables in "G()", "time_step()" and "calculate_forcing()"             
clear all; 

dbstop if error
format short g

addpath('MEX');
addpath(genpath('utils'));
addpath(genpath('plotting'));

% % FIXME: diary behavior is such that if renamed it's still active diary!!
diary off; diary off; diary on; diary off; diary on; diary on;

% Clean up for threads, more complex and much slower than expected.

timer_total = tic;

% killAndCleanThreads();  % Leftover threads and dumps cause BIG trouble
% toc(timer_total)

%%%%%%

% Setup the big picture parts of a simulation and/or NK solution.
% FIXME: Someday, when we know what inputs need to be, put all this a file

tName = tracer_names(0);    % no CISO tracers
selection = [ ...
    %     find(strcmp(tName,'SiO3')),...
    %     find(strcmp(tName,'DIC')),...
    find( strcmp(tName,'O2') ) ];
forwardIntegrationOnly = 0;
ck_years = 1;   % Newton-Kryrlov -requires- 1 year, but might want to run long time_step_hr = 3;
time_step_hr = 3;
logTracers = 0;
yearsBetweenRestartFiles = 1;
captureAllSelectedTracers = 0;

% DEBUG stuff
% logTracers = 1;
% ck_years = 1;
% time_step_hr = 12; % FAST debug
% yearsBetweenRestartFiles = 1;
% % time_step_hr = 6912/60/60;

% captureAllSelectedTracers = 1;

%%%%%%

marbl_file = 'Data/marbl_in'; % MARBL chemistry and other constants.

%%%%%% INput restart file

% start_yr = 0; inputRestartFile = 'Data/passive_restart_init.mat'; % from netCDF 5/25/22
start_yr =  70; inputRestartFile = 'Data_GP/restart_70_integrate_from_0.mat';
start_yr = 260; inputRestartFile = 'Data_GP/restart_260_integrate_from_0.mat';
% start_yr = 4101; inputRestartFile = 'Data/InputFromAnn/restart4101.mat';

fprintf('%s.m: Reading OFFline input restart file with tracers and transports: %s\n', mfilename, inputRestartFile);
load(inputRestartFile,'sim','MTM');

% We just over wrote sim struct, but now we can save the file stuff in it.
% Need input restart name to make our output...

sim.forwardIntegrationOnly  = forwardIntegrationOnly ;
sim.inputRestartFile        = inputRestartFile; 
sim.start_yr                = start_yr;
sim.selection               = selection;
sim.captureAllSelectedTracers = captureAllSelectedTracers;

sim.logTracers = logTracers;
sim.logDiags   = and (0, sim.logTracers) ; % Usually no diags..

% FIXME: lots of old code floating around...
clear inputRestartFile start_yr selection captureAllSelectedTracers logTracers 

sim.checkNeg = 0;

%%%%%% OUTput restart file

% "restart file" -FROM- this OFFline sim, if for restart on Mac, not CESM.
% Another hack is needed to move results from NK here back to CESM.

sim.outputRestartDir = myRestartDir(ck_years);
disp(['Results will be saved in directory ', sim.outputRestartDir]); disp (' ');
[status, msg, msgID] = mkdir(sim.outputRestartDir); 
if status ~=1
    disp(msg); disp(msgID); disp(' ')
    keyboard
else
    clear status msg msgID
end


sim.yearsBetweenRestartFiles = yearsBetweenRestartFiles;       % YEARS between drops of restart

%%%%%%

% parallel is hard to debug, but 2x faster

sim.runInParallel = 1;
if (sim.runInParallel)
    sim.number_of_threads = 4; end

%%%%%%

% In past I debuged MARBL Carbon isotopes. "lciso_on", and that stuff, it
% probably still works but they makes everything bigger and mch slower.

sim.lciso_on = 0;   % run with Carbon Isotopes ??

%%%%%%

sim.epsilon = -sqrt(eps);

%%%%%%

% I think it's impossible to debug anything with out looking at time step
% data. If you can debug without intermediate results just turn off
% logging.
% %
% However MARBL diagnostics, while occasionally invaluable are HUGE and
% slow the sim down by at least a factor 2x, require huge memry etc; so I
% very rarely capture them...

if (or (sim.logDiags, sim.logTracers))
    disp('Setting up to log tracers and possibly MARBL diags...')
    % need grid dimensions for many things. Starting with possible log files...
    [~, sim.domain.wet_loc, sim.domain.iwet_FP, ~, sim.domain.bottom_lvl, ~]...
        = createIwetEtc(sim.domain.M3d);

    % Possibly Record everything at a single location, water level. But which?
    % FIXME: need a clean way to convert (lat,lon) to (iLat,iLon), but in the
    % mean time; iterate thise by hand...
    %
    % iLat = 50; iLon = 81; iLvl = 10;    % Galapagos(0.3N, 108.7W)   iFp = 6496
    % iLat = 77; iLon = 60; iLvl = 10;    % NorPac  (45.7N, 176.8E)   iFp = 4629
    % iLat = 50; iLon = 61; iLvl = 1;     % Dateline  ( 0.3N,  179.3E)   iFp = 4702
    % iLat = 49; iLon = 11; iLvl = 10;    % Zulu =  ( 0.3N, 0.7E)     iFp = 1049
    % iLat = 58; iLon = 50; iLvl = 10;    % Palau = ( 5.6N, 139.7E)   iFp = 3704
    % iLat = 50; iLon = 28; iLvl = 10;    % IO      ( 0.3N,  50.5E)   iFp = 2080
%     iLat = 57; iLon =  3; iLvl = 10;    % AF 447 =  ( 4.7N, -29.5E)     iFp = 1049
%     iLat = 20; iLon =  95; iLvl = 4;    % "-48" =  ( -45.695N, -58.3E)     iFp = 31045 iCol 7462
    iLat = 2; iLon =  95; iLvl = 1;    % 7445
    % Check that! Make a map!
    % first get iFp on level 1, Simpy put: on level 1, iFp = iCol...

    iFp = coordTransform_xyz2fp(iLat, iLon, 1, sim);
    [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(iFp, sim, 999); title('Time Series Localtion')
%     [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(iFp, sim); 
    sim.time_series_loc = iFp ;
    % % % ... now set the level
    sim.time_series_lvl = iLvl;
    disp(['Time series(loc,lvl) = (', num2str(sim.time_series_loc), ', ', num2str(sim.time_series_lvl),')']);
    [~,~,~, lat, lon, ~] = col2latlon(sim, sim.time_series_loc);
    disp(['Time series(lat,lon, depth) = (', num2str(lat,'%.1f'), ' N, ', num2str(lon,'%.1f'),' E, ',num2str(sim.grd.zt(sim.time_series_lvl),'%.1f'),' m))'])
    disp(' ')
    clear iLat iLon iLvl lat lon iFp
else
    % avoid messy code in parallel;
    % just set a default legal array idx
    sim.time_series_loc = 1 ;
    sim.time_series_lvl = 1;
end

%%%%%% End of "inputs"
toc(timer_total)

% Stuff below is not a simulation input. It is code to setup grids, etc a
% forward integration, or NK(), or...
%
[sim, bgc, ~, time_series, forcing] = init_sim(marbl_file, sim.inputRestartFile, sim, ck_years, time_step_hr);
if (or (sim.logDiags, sim.logTracers))
    getMemSize(time_series,1e3);
end

% DEBUG: run on a few locations, rather global, MUCH faster debug something
% sim.domain.num_wet_loc = 1; % comment out too run entire world

% scale units to unity, using values of global mean of initial
% need a column vector

sim = calc_global_moles_and_means(bgc, sim);

clear forwardIntegrationOnly marbl_file
%%%%%%

toc(timer_total)

% =============== This is the NK solver code ================

fprintf('\n%s.m: Start Newton-Krylov solver: %s\n',mfilename,datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
tic;

% NK always using 1 year integration

if ck_years ~= 1
    disp('ERROR: NK requires ck_years = 1')
    keyboard
end

% bgc of default

% c0 = bgc2nsoli(sim, bgc.tracer);

% ALWAYS punt "alt" methods that do NOT depend on any other tracers...

idx = ismember(sim.selection, [9,11]); 
sim.selection(idx)=[];

sim.selection = unique(sort(sim.selection));
if (min(sim.selection)<1) || (max(sim.selection)>32)
    keyboard   % bogus tracer selected
end
cstr = tName(sim.selection)'; fprintf(1,'Selected tracers: '); fprintf(1,'%s ',cstr{:}); fprintf('\n');

%%%%
% keyboard
disp('Parameters nsoli()...')

nsoli.lmeth  = 2;               % method 2 = GMRES(m)
nsoli.atol   = 1e-3;
nsoli.rtol   = 1e-9;            % stop when norm is less than atol+rtol*norm of init_resid as seen by nsoli
nsoli.tol    = [nsoli.atol,nsoli.rtol];     % [absolute error, relative tol]
nsoli.etamax = 0.9;             % maximum error tol for residual in inner iteration, default = 0.9
nsoli.maxit  = 30;              % maximum number of nonlinear iterations (Newton steps) default = 40
nsoli.maxitl = 10;              % maximum number of inner iterations before restart in GMRES(m), default = 40;
                                % also number of directional derivative calls, also num of gmres calls
nsoli.restart_limit = 10;       % max number of restarts for GMRES if lmeth = 2, default = 20;
parms  = [nsoli.maxit, nsoli.maxitl, nsoli.etamax, nsoli.lmeth, nsoli.restart_limit];


fprintf('%s.m: Starting nsoli()...\n',mfilename);
timer_loop = tic;

tot_t = sim.dt*sim.num_time_steps;  % used only for debug output
disp([num2str(sim.num_time_steps),' time steps, each time step is ', num2str(sim.dt/sim.const.sec_h), ...
    ' (h), simulating ', num2str(tot_t /sim.const.sec_y, '%1.3f'),' years'])

% DEBUG: add a known spike to see if it gets removed by Nsoli()spike_size = 1234.567;
spike_size = 0;
% spike_size = 1234.567;
fprintf('\n**** %s.m: Adding a spike of %g to water parcel (%d,%d,%d) ****\n\n', mfilename, spike_size, sim.time_series_loc, sim.time_series_lvl,sim.selection);
bgc.tracer(sim.time_series_loc,sim.time_series_lvl,sim.selection) = spike_size +bgc.tracer(sim.time_series_loc,sim.time_series_lvl,sim.selection);


c0 = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
sz = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
c = reshape(c0, sz);    

% use nsoli() only on selected tracers

x0 = c(:,sim.selection);    % initial condition for Nsoli()
x0 = x0(:);             % unitless

% FIXME: use average of something or instanteous or ...
% average_tendency = load('avgTend.mat','average_tendency').average_tendency;
% average_tracer   = load('avgTracer.mat','average_tracer').average_tracer;
% bgc_0 = average_tracer;  % [272000,1]
% 
% Calculate Jacobian of average tracer tendency; Not a slow as it sounds...
%
% Note: MARBL is not exactly linear in the tracer, but not so far off from
% that. So while the average of the tendency is -NOT- the tendency of the
% average tracer, its not far off. Both the average of the tracer and the
% average tendency are captured so that the error from this assumption can
% be determined.
% 
% % % Calculate J or load saved?
% % if (0)
% %     J = calc_J_full(@calc_f, packMarbl(bgc.tracer, sim.domain.iwet_JJ), sim, bgc, time_series);
% %     Q_inv = calc_Q_inv(J, bgc, sim);
% %     save('/Users/jj/Desktop/UCI/MARBL/MARBL_v4/J.mat','J');
% %     save('/Users/jj/Desktop/UCI/MARBL/MARBL_v4/Q_inv.mat','Q_inv');
% %     %     keyboard
% % else
% %     load('/Users/jj/Desktop/UCI/MARBL/MARBL_v4/J.mat','J');
% % %     load('/Users/jj/Desktop/UCI/MARBL/MARBL_v4/Q_inv.mat','Q_inv');
% % end
% 
% Q = cycle period *fPrime
%
% fPrime = Transport +J
% fPrime = blkdiag(tmp_TR{:}) +J;
% getMemSize(fPrime,1e3);
%   but don't actually need fPrime, just fPrime*T
%
% Calculate inverse of Q for each water column: all levels all tracers.
% keyboard
% 
% % if (0)
% %     Q_inv = calc_Q_inv(J, bgc, sim);
% %     save('/Users/jj/Desktop/UCI/MARBL/MARBL_v4/Q_inv.mat','Q_inv');
% % else
% %     load('/Users/jj/Desktop/UCI/MARBL/MARBL_v4/Q_inv.mat','Q_inv');
% % end
% 
% P = speye(size(Q)) -Q; P is never used. Do NOT use mfactor() to solve Q
% 
% PQ_inv = Q_inv -speye(size(Q_inv));
% 
% getMemSize(J,1e3);
% getMemSize(Q_inv,1e3);
% getMemSize(PQ_inv,1e3);
% clear J Q_inv
% keyboard
% plot_log(128, "Global mean", myDepth, sim.globalMean_lvl, nameUnits(idx), idx, true);
disp('FIXME: Set PQ_inv = 1 for now, but we need preconditioner...')
PQ_inv = 1;

[sol,it_hist,ierr,x_hist] = Cnsoli(x0, @(x0,Npt) calc_G(x0,Npt,c0,sim,bgc,time_series,forcing,MTM,PQ_inv), nsoli.tol, parms);

% convert solution to include all tracers, even ones not optimized x = reshape(c0,sz);

x0 = reshape(c0, sz);
x0(:,sim.selection) = reshape(sol,[sz(1),numel(sim.selection)]);
x0 = x0(:);
bgc.tracer(:) = nsoli2bgc(sim, bgc, x0);


disp('Nsoli() finished...')
keyboard
elapsedTime_all_loc = toc(timer_loop);
disp(' ');
disp(['Runtime: ', num2str(elapsedTime_all_loc, '%1.0f'),' (s) or ', num2str(elapsedTime_all_loc/60, '%1.1f'), ' (m)'])
disp(['Runtime per location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps/sim.domain.num_wet_loc*1000, '%1.2f'), ' (ms) MARBL, advection, diffusion, mfactor()'])
disp(['Runtime all location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps, '%1.2f'),                    ' (s)  MARBL, advection, diffusion, mfactor()'])
disp(['Runtime all location per sim year : ', num2str(elapsedTime_all_loc/60/1440/tot_t*sim.const.sec_y, '%1.2f'), ' (d/y_sim)'])
disp(['Simulation speed: ', num2str(tot_t/elapsedTime_all_loc/sim.const.days_y, '%1.1f'), ' (sim y/d) aka (SYPD)'])

keyboard
% FIXME: need to save workspace?!
save_timer = tic; disp('Saving (possibly) large workspace file...'); save(strcat(sim.data_dir+'/Logs/last_run.mat')); toc(save_timer);

disp(['Log file of one location for all time steps uses ',num2str(getMemSize(time_series)/1024/1024, '%1.1f'),' MB'])
disp(['Log file of one location uses ',num2str(getMemSize(time_series)/1024/sim.num_time_steps, '%1.1f'),' KB per time step'])




% %%% FIXME: Hackery for globals needs to go aways
% 
% global sim bgc time_series forcing MTM sim.selection
% global PQ_inv
% 
%%%%%%
%
% Cleanup left over threads, old variables etc, etc, etc
%
% MATLAB "analyze code" says "clear all" or "mex" is slow, but to avoid
% strange hangs and crashes if MEX does not end gracefully on previous run,
% you have to clear mex. "clear all" also removes any functions that have
% been run and therefore parsed etc, but who cares?
% clear mex; clear variables; clear global;

