% "main" of cyclostationary transport version of MARBL,
%
% Try to solve MARBL using "Newton Like" methods from Kelley
%
%        https://ctk.math.ncsu.edu/newtony.html
%
%   matlab -nodisplay -nodesktop -nosplash -noFigureWindows -logfile batch.txt < marbl_nsoli.m &
%
% MARBL code being used is https://github.com/marbl-ecosys/MARBL
%
% function marbl()
% "main" of cyclostationary transport version of MARBL,

fprintf('%s.m: Start at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
% Need this to clear the "persistent" variables in "G()", "time_step()" and "calculate_forcing()"
clear all;
dbstop if error
format short g
addpath('MEX');
addpath(genpath('utils'));
addpath(genpath('plotting'));
% dname = sprintf('%s/../',myDataDir());
% addpath(dname);
% % FIXME: diary behavior is such that if renamed it's still active diary!!
diary off; diary off; diary on; diary off; diary on; diary on;
% Clean up for threads, more complex and much slower than expected.
% killAndCleanThreads();  % Leftover threads and dumps cause BIG trouble
timer_total = tic;

%%%%%%

% Setup big picture parts of a simulation and/or NK solution.

% always need a selected tracer! For plot time series, or solve
% Most of these work very well in the single tracer solution...
tName = tracer_names(0);    % no CISO tracers
% tracer_str = 'DOCr';
tracer_str = 'DOP';
% tracer_str = 'ALK';
selection = [ find( strcmp(tName,tracer_str) ) ];
% selection = [ ...
%     find( strcmp(tName,'SiO3') ) ];     % #3
%     find( strcmp(tName,'DOCr') ) ];     % #17

forwardIntegrationOnly = 0; % 1 -> no NK just fwd integration
recalculate_PQ_inv     = 1; % recalculate J, PQ,PQ_inv or load file
% remember that "sol" is an x0 value. First relax iteration gives same x1
% as sol run. to be useful num_relax_iterations >= 2
num_relax_iterations   = 2; % 0 means no relax steps, just use NK sol
num_forward_years      = 2; % if fwd only, this is num fwd;
% else this is num fwd after relax , after sol
logTracers                = 0;
captureAllSelectedTracers = 0;
yearsBetweenRestartFiles  = 1;

% FIXME: Someday, when we know what inputs need to be, put all this a file
time_step_hr = 3;
phi_years    = 1;      % NK always using 1 year integration

% DEBUG stuff
logTracers = 1;
time_step_hr = 12; % FAST debug

%%%%%%
marbl_file = 'Data/marbl_in'; % MARBL chemistry and other constants.
%%%%%%

% Input restart file

% start_yr = 4101;inputRestartFileStem = 'Data/InputFromAnn/restart4101.mat';
% start_yr =   0; inputRestartFileStem = 'Data/passive_restart_init.mat'; % from netCDF 5/25/22
% start_yr = 260; inputRestartFile = 'Data_GP/restart_260_integrate_from_0.mat';
% start_yr = 2525; inputRestartFileStem = 'restart_0_1_output/restart_260_NH4_x0_sol.mat';
start_yr = 12345; inputRestartFileStem = 'restart_0_1_output/restart_3536_Fe_x1.mat';

inputRestartFile = strcat(myDataDir(), inputRestartFileStem);
fprintf('%s.m: Reading OFFline input restart file with tracers and transports: %s\n', mfilename, inputRestartFile);
% load() does NOT need file extension, but copy() does. sigh
if ~isfile(inputRestartFile)
    error("missing or typo in name of inputRestartFile")
    keyboard
end
load(inputRestartFile,'sim','MTM');

% We just over wrote sim struct, so now we can save stuff in it again. sigh

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

% Need input restart name to make our output ???
%
% "restart file" -FROM- this OFFline sim for restart on Mac, not CESM.
% Another hack is needed to move results from NK here back to CESM.

sim.yearsBetweenRestartFiles = yearsBetweenRestartFiles;       % YEARS between drops of restart

sim.outputRestartDir = strcat(myDataDir(),'restart_0_1_output/');
disp(['Results will be saved in directory ', sim.outputRestartDir]); disp (' ');
[status, msg, msgID] = mkdir(sim.outputRestartDir);
if status ~=1
    disp(msg); disp(msgID); disp(' ')
    keyboard
else
    clear status msg msgID
end

%%%%%%

% parallel is hard to debug, but 2x faster

sim.runInParallel = 1;
if (sim.runInParallel)
    sim.number_of_threads = 4; % only 4 on laptop or 12 on GP supported
end

%%%%%%

% In past I debuged MARBL Carbon isotopes. "lciso_on", and that stuff, it
% probably still works but they makes everything bigger and mch slower.

sim.lciso_on = 0;   % run with Carbon Isotopes ??

%%%%%%

sim.epsilon = -sqrt(eps);

%%%%%%

sim = setPeek(sim);

%%%%%% End of "inputs"

% Stuff below is not a simulation input. It is code to setup grids, etc a
% forward integration, or NK(), or...
%

[sim, bgc, ~, time_series, forcing] = init_sim(marbl_file, sim.inputRestartFile, sim, phi_years, time_step_hr);

tot_t = sim.dt*sim.num_time_steps;  % used only for debug output
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

% =============== This is the NK solver code ================

toc(timer_total)

fprintf('\n%s.m: Start Newton (Broyden Method) solver: %s\n',mfilename,datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
timer_PQ_init_solve_relax_fwd = tic;


% NK always using 1 year integration

if phi_years ~= 1
    disp('ERROR: NK requires phi_years = 1')
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
cstr = tName(sim.selection)';
fprintf('%s: Selected tracers: %s \n', mfilename, cstr{:});

%%%%
% % DEBUG: add a known spike to see if it gets removed by Nsoli()spike_size = 1234.567;
% spike_size = 0;
% % spike_size = 1234.567;
% fprintf('\n**** %s.m: Adding a spike of %g to water parcel (%d,%d,%d) ****\n\n', mfilename, spike_size, sim.time_series_loc, sim.time_series_lvl,sim.selection);
% bgc.tracer(sim.time_series_loc,sim.time_series_lvl,sim.selection) = spike_size +bgc.tracer(sim.time_series_loc,sim.time_series_lvl,sim.selection);
%%%%

c0 = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
sz = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
c = reshape(c0, sz);

% Solve only on selected tracers

x0 = c(:,sim.selection);    % initial condition for Nsoli()
x0 = x0(:);                 % unitless

if( sim.forwardIntegrationOnly )
    fprintf('%s.m: forward integration ONLY\n',mfilename);
else
    % Solve for selected tracer
    if(recalculate_PQ_inv)
        PQ_inv = calc_PQ_inv(sim, bgc, time_series, forcing, MTM);
    else
        fprintf('\n%s.m: Loading ~30 GB(!) mfactored preconditioner PQ_inv from %s...\n', mfilename, strcat('sol/',string(tName(sim.selection)),'_QJ'))
        tStart = tic;
        %     keyboard
        load (strcat(myDataDir(),'sol/',strjoin(tName(sim.selection)),'_QJ'), 'PQ_inv')

        fprintf('%s.m: %1.0f (s) to init sim and load PQinv \n',mfilename, toc(tStart));
    end % calculate or load PQ_inv

    fprintf('%s: Parameters nsoli()... \n', mfilename)
    maxit  = 7;            % maximum number of nonlinear iterations (Newton steps) default = 40
    maxitl = 2;            % maximum number of Broyden iterations before restart, so maxdim-1 vectors are stored default = 40

    % used only by nsoli()
    etamax = 0.9;           % maximum error tol for residual in inner iteration, default = 0.9
    lmeth  = 2;             % Nsoli() method 2 = GMRES(m), not used by brsola().
    restart_limit = 10;     % max number of restarts for GMRES if lmeth = 2, default = 20;

    % first 2 of these parms are used by brsola, reat are specific to nsoli
    parms  = [maxit,maxitl,  etamax,lmeth,restart_limit];

    % stop when norm is less than atol+rtol*norm of init_resid as seen by nsoli or brsola
    %
    % atol   = sqrt(eps);     % sum of the squares in (1/s), IOW average error = 1/sec_y/379,913 ~1e-11 years
    % rtol   = 1e-3;          % stop when norm is less than atol+rtol*norm of init_resid as seen by nsoli
    %
    atol   = sqrt(eps);     % sum of the squares in (1/s), IOW average error = 1/sec_y/379,913 ~1e-11 years
% rather than get accurate solution, loop over tracers and try to reduce G
% to 1% of starting value, while looping over all the tracers, and then
% repeat until to get final very accurate result where G is sqrt(eps)
    [r,G, x1] = calc_G(x0,c0,sim,bgc_sol,time_series,forcing,MTM,PQ_inv);
    rtol   = G *1e-2; % stop if norm(drift) < 1ppm of G(x0)

    tol    = [atol,rtol];   % [absolute error, relative tol]
    [sol,it_hist,ierr,x_hist] = brsola(x0, @(x) calc_G(x,c0,sim,bgc,time_series,forcing,MTM,PQ_inv), tol, parms);

    % remember that "sol" is an x0 value...
    %
    % sol = [379913]
    % c0  = [379913, 32]
    % tracer = [7881, 60, 32]
    nsoli_sol = replaceSelectedTracers(sim, c0, sol, sim.selection);    % [12157216] = [379913*32]
    bgc_sol = bgc;
    bgc_sol.tracer = nsoli2bgc(sim, bgc_sol, nsoli_sol);

    myRestartFile = sprintf('%s/restart_%s_sol.mat', sim.outputRestartDir, strjoin(tName(sim.selection)));
    [sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile);


    sol_fname = sprintf('%s/sol_%s_ierr_%d', sim.outputRestartDir, string(tName(sim.selection)), ierr);
    fprintf('%s.m: Saving sol, ierr, it_hist, x_hist in %s"...\n', mfilename, sol_fname);
    save(sol_fname, 'sol', 'ierr', 'it_hist', 'x_hist', '-v7.3','-nocompression');

    fprintf('\n%s.m: Now we "relax" or forward integtration single variable solution a few iterations: %s\n',mfilename,datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
    fprintf('%s.m: %d time steps, each time step is %1.1f (h), simulating %1.1f years\n', ...
        mfilename, sim.T/sim.dt, sim.dt/sim.const.sec_h, tot_t /sim.const.sec_y)

    % remember that "sol" is an x0 value. First relax iteration gives same
    % x1 of sol. to be useful num_relax_iterations >= 2
    x = sol;
    %     x_histx = x;
    %     r_hist = zeros(size(x));
    %     it_histx = zeros(num_relax_iterations,1);
    for itc = 1:num_relax_iterations

        fprintf("\n%s.m: starting relaxation year #%d of %d\n", mfilename, itc, num_relax_iterations)

        % Note x = x1, not "x0" which is normal thing for sol iterations
        [r,G, x1] = calc_G(x,c0,sim,bgc_sol,time_series,forcing,MTM,PQ_inv);

        % DEBUG
        fprintf('%s.m: itc %d norm(G) = %g\n',mfilename, itc, norm(G))
        fprintf('%s.m: itc %d norm(r) = %g\n',mfilename, itc, norm(r))
        %         x_histx = [x_histx,x];
        %         r_hist  = [r_hist, r];
        %         it_histx(itc,1) = norm(r);

        x = x1;

    end % relax loop

    % If we relaxed sol, x is relaxed x1.
    %   if we did NOT relax, then x = sol = x0

    nsoli_relax = replaceSelectedTracers(sim, c0, x, sim.selection);
    bgc_relax = bgc;
    bgc_relax.tracer = nsoli2bgc(sim, bgc_relax, nsoli_relax);

    bgc = bgc_relax;    % final answer or input to fwd integration...

    myRestartFile = sprintf('%s/restart_%s_relax_x1.mat', sim.outputRestartDir, strjoin(tName(sim.selection)));
    [sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile);
end

% Next! allow ALL tracers, not just selection, to "relax' to solution.
%    do pure forward integration for a while...
years_gone_by = 0;
bgc_fwd = bgc;
for fwd_itc = 1:num_forward_years
    fprintf("\n%s.m: starting forward integrate year #%d of %d\n", mfilename, fwd_itc, num_forward_years)
    [sim, bgc_fwd, time_series, tracer_0] = phi(sim, bgc_fwd, time_series, forcing, MTM);

    sim.start_yr = sim.start_yr+1;

    years_gone_by = years_gone_by +1;
    if mod(years_gone_by, sim.yearsBetweenRestartFiles) == 0    % This runs after last time step of every 10 y
        % years gone by is actually the number of years-1 spent in phi.
        myRestartFile = sprintf('%s/restart_%s_fwd_x1.mat', sim.outputRestartDir, strjoin(tName(sim.selection)));
        [sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile);
    end
end % fwd loop
bgc = bgc_fwd; % this is in fact my final answer!

elapsedTime_all_loc = toc(timer_PQ_init_solve_relax_fwd);
disp([mfilename,' finished...'])
disp(' ');
disp(['Runtime: ', num2str(elapsedTime_all_loc, '%1.0f'),' (s) or ', num2str(elapsedTime_all_loc/60, '%1.1f'), ' (m)'])
disp(['Runtime per location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps/sim.domain.num_wet_loc*1000, '%1.2f'), ' (ms) MARBL, advection, diffusion, mfactor()'])
disp(['Runtime all location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps, '%1.2f'),                    ' (s)  MARBL, advection, diffusion, mfactor()'])
disp(['Runtime all location per sim year : ', num2str(elapsedTime_all_loc/60/1440/tot_t*sim.const.sec_y, '%1.2f'), ' (d/y_sim)'])
disp(['Simulation speed: ', num2str(tot_t/elapsedTime_all_loc/sim.const.days_y, '%1.1f'), ' (sim y/d) aka (SYPD)'])

% keyboard
% years gone by is actually the number of years-1 spent in phi.
[sim, bgc] = saveRestartFiles(sim, bgc, bgc.tracer, myRestartFile);
% FIXME: need to save workspace?!
logDir = strcat(sim.outputRestartDir,'/Logs/');
if ~exist(logDir, 'dir')
    mkdir(logDir)
end
save_timer = tic; disp('Saving (possibly) large workspace file...'); save(strcat(logDir,'last_run.mat'),'-v7.3','-nocompression'); toc(save_timer);

disp(['Log file of one location for all time steps uses ',num2str(getMemSize(time_series)/1024/1024, '%1.1f'),' MB'])
disp(['Log file of one location uses ',num2str(getMemSize(time_series)/1024/sim.num_time_steps, '%1.1f'),' KB per time step'])
