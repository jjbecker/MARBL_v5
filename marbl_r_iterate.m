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
logTracers = 0;
ck_years = 1;
time_step_hr = 12; % FAST debug
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
    iLat = 20; iLon =  95; iLvl = 4;    % "-48" =  ( -45.695N, -58.3E)     iFp = 31045 iCol 7462
    %     iLat = 2; iLon =  95; iLvl = 1;    % 7445
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

% nsoliParam.lmeth  = 2;               % method 2 = GMRES(m)
% nsoliParam.atol   = 1e-3;
% nsoliParam.rtol   = 1e-9;            % stop when norm is less than atol+rtol*norm of init_resid as seen by nsoli
% nsoliParam.tol    = [nsoliParam.atol,nsoliParam.rtol];     % [absolute error, relative tol]
% nsoliParam.etamax = 0.9;             % maximum error tol for residual in inner iteration, default = 0.9
% nsoliParam.maxit  = 30;              % maximum number of nonlinear iterations (Newton steps) default = 40
% nsoliParam.maxitl = 10;              % maximum number of inner iterations before restart in GMRES(m), default = 40;
%                                 % also number of directional derivative calls, also num of gmres calls
% nsoliParam.restart_limit = 10;       % max number of restarts for GMRES if lmeth = 2, default = 20;
% parms  = [nsoliParam.maxit, nsoliParam.maxitl, nsoliParam.etamax, nsoliParam.lmeth, nsoliParam.restart_limit];


fprintf('%s.m: Starting r iteration for tracer selection = [%d]...\n',mfilename, sim.selection);
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
x0 = x0(:);                 % unitless

disp('FIXME: Set PQ_inv = 1 for now, but we need preconditioner...')
% %     J = calc_J_full(@calc_f, x0, sim, bgc, time_series);
PQ_inv = 1;


% From NK-EXAMPLES: make a sparse operator that restores the surface to zero with a time scale of tau
msk  = sim.domain.M3d;     % wet == 1, dry == 0
iwet = sim.domain.iwet_FP;
% tau = 24 * 60^2;                % (sec/d)
% tau = 7*sim.const.sec_d;        % (sec/wk)
tau = sim.const.sec_h;          % (sec/hr)
temp = msk;
temp(:,:,2:end) = 0;
R   =  d0( temp(iwet) / tau );   % (1/sec)

Q =  MTM(1).A + MTM(1).H + MTM(1).D - R;
for k = 2:12;
    Q = Q + MTM(k).A + MTM(k).H + MTM(k).D - R;
end
%         T = 12*num_step_per_month*dt;
Q = Q/12; % annually averaged transport and surface restoring (aka birth)
fprintf('  Factoring the preconditioner...');
tic;
FQ = mfactor(sim.T*Q);
%         sim.FQ = FQ;
%         toc
%         x0_age = sim.domain.M3d + nan;
%         x0_age(sim.domain.iwet_FP) = -mfactor(FQ,0*iwet+sim.T);
%
% %     J = calc_J_full(@calc_f, x0, sim, bgc, time_series);
PQ_inv = FQ;
toc


% for jj=1:5
%     J_wrong = calc_J_simplified(@calc_f, x0, sim, bgc, time_series, forcing);
% %     J_o2 = calc_J_O2(@calc_f, x0, sim, bgc, time_series, forcing);
% end
num_r_iterations = 40;
x = x0;
% x = load('Data/restart_0_1_output_rIter/G_30.mat', 'x0').x0;
% load('Data_GP/saveRes_25.mat', 'x_hist');
% x = x_hist(:,25);
clear x_hist

x_hist = x;
r_hist = zeros(size(x));
it_histx = zeros(num_r_iterations,1);
Npt = -1;

% num_r_iterations = 0;
for itc = 1:num_r_iterations
    % FIXME: note "x" not "x0"

    [r,G] = calc_G(x,c0,sim,bgc,time_series,forcing,MTM,PQ_inv);
% disp('FIXME: hacking r = -G')
% r = -G;
    fnrm = norm(r);
    fprintf('%s.m: itc %d norm(G) = %g\n',mfilename, itc, norm(G))
    fprintf('%s.m: itc %d norm(r) = %g\n',mfilename, itc, fnrm)

    x_hist = [x_hist,x];
    r_hist = [r_hist,r];
    it_histx(itc,1) = norm(r);

    % recursion:    x <- x -r;
    % w = 1;    % x = 2x -phi if iterating with G
    w = 0.9;
    x = x -w*r;     % e.g. x = x +w*f(x) - w*x = (1-w)x + wf(x)

end
sol = x;

% % Convert the solution into a standard MARBL format
% %
% % [sol, it_hist, ierr, x_hist] = Cnsoli(x0, @(x0) calc_G(x0,c0,sim,bgc,time_series,forcing,MTM,PQ_inv), nsoliParam.tol, parms);
% [sol, it_hist, ierr, x_hist] = nsoli(x0, @(x0) calc_G(x0,c0,sim,bgc,time_series,forcing,MTM,PQ_inv), nsoliParam.tol, parms);
% 
% x0_bgc = replaceSelectedTracers(sim, c0, sol, sim.selection);
% bgc.tracer = nsoli2bgc(sim, bgc, x0_bgc);   % marbl format
%


% Run the solution forward a few years to let the non-optimized tracers
% "relax" to the new values.

for num_r_relax_iterations = 1:10
    fprintf("%s.m: starting relaxation year #%d\n", mfilename, num_r_relax_iterations)
    [sim, bgc, time_series] = phi(sim, bgc, time_series, forcing, MTM);
end

disp('r_iterate() finished...')

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
