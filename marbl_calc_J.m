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

addpath '/Users/jj/Desktop/UCI/MARBL/MARBL_v5/sol'

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
% selection = [ ...
%     find( strcmp(tName,'SiO3') ) ];     % #3
%     find( strcmp(tName,'DOCr') ) ];     % #17
%
% Most of these work very well in the single tracer solution...
tracer_str = 'diazC';

selection = [ find( strcmp(tName,tracer_str) ) ];
forwardIntegrationOnly = 0;
ck_years = 1;   % Newton-Kryrlov -requires- 1 year, but might want to run long time_step_hr = 3;
time_step_hr = 3;
logTracers = 0;
yearsBetweenRestartFiles = 1;
captureAllSelectedTracers = 0;

% DEBUG stuff
logTracers = 1;
ck_years = 1;
time_step_hr = 12; % FAST debug
% captureAllSelectedTracers = 1;

%%%%%%

marbl_file = 'Data/marbl_in'; % MARBL chemistry and other constants.

%%%%%% INput restart file

% start_yr = 4101;inputRestartFile = 'Data/InputFromAnn/restart4101.mat';
% start_yr =   0; inputRestartFile = 'Data/passive_restart_init.mat'; % from netCDF 5/25/22
% start_yr = 260; inputRestartFile = 'Data_GP/restart_260_integrate_from_0.mat';
  start_yr = 260; inputRestartFile = 'Data/restart_0_1_output/restart_260_NH4_x0_sol.mat';
%   start_yr = 260; inputRestartFile = 'Data/restart_0_1_output/restart_260_DiazC_x0.mat';

fprintf('%s.m: Reading OFFline input restart file with tracers and transports: %s\n', mfilename, inputRestartFile);
% load() does NOT need file extension, but copy() does. sigh
if ~isfile(inputRestartFile)
    error("missing or typo in name of inputRestartFile")
    keyboard
end
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

sim.outputRestartDir = myDataDir(ck_years);
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

sim = setPeek(sim);

%%%%%% End of "inputs"
toc(timer_total)

% Stuff below is not a simulation input. It is code to setup grids, etc a
% forward integration, or NK(), or...
%

[sim, bgc, ~, time_series, forcing] = init_sim(marbl_file, sim.inputRestartFile, sim, ck_years, time_step_hr);

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

toc(timer_total)

% =============== This is the NK solver code ================

fprintf('\n%s.m: Start Newton (Broyden Method) solver: %s\n',mfilename,datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
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


selection = sim.selection;

c0 = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
sz = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
c = reshape(c0, sz);

if(0)
    [PQ_inv, J_FP] = calc_PQ_inv(sim, bgc, time_series, forcing, MTM);
else
    fprintf('\n%s.m: Loading ~30 GB(!) mfactored preconditioner PQ_inv from %s...\n', mfilename, strcat('sol/',string(tName(sim.selection)),'_QJ'))
    tStart = tic;
%     keyboard
    load (strcat(string(tName(sim.selection)),'_QJ'), 'PQ_inv')

    elapsedTime = toc(tStart);
    fprintf('%s.m: %1.0f (s) to init sim and load PQinv \n',mfilename, toc(tStart));
end

fprintf('%s: Parameters nsoli()... \n', mfilename)

% stop when norm is less than atol+rtol*norm of init_resid as seen by nsoli or brsola
% 
% atol   = eps;           % sum of the squares in (1/s), IOW average error = 1/sec_y/379,913 = 8e-14 years
% rtol   = 1e-7;          % stop when norm is less than atol+rtol*norm of init_resid as seen by nsoli
% 
% atol   = 1;
% rtol   = 1e-7;          % stop when norm is less than atol+rtol*norm of init_resid as seen by nsoli
% 
rtol   = 0;             % unfortunatley this is used with the initial drift, which is completely unknown                                
atol   = norm(x0)*1e-6; % stop if norm(drift) < 1ppm of norm(x0)

tol    = [atol,rtol];   % [absolute error, relative tol]
maxit  = 40;            % maximum number of nonlinear iterations (Newton steps) default = 40
maxitl = 15;            % maximum number of Broyden iterations before restart, so maxdim-1 vectors are stored default = 40

% used only by nsoli()
etamax = 0.9;           % maximum error tol for residual in inner iteration, default = 0.9
lmeth  = 2;             % Nsoli() method 2 = GMRES(m), not used by brsola().
restart_limit = 10;     % max number of restarts for GMRES if lmeth = 2, default = 20;

% first 2 of these parms are used by brsola, reat are specific to nsoli
parms  = [maxit,maxitl,  etamax,lmeth,restart_limit]; 

[sol,it_hist,ierr,x_hist] = brsola(x0, @(x) calc_G(x,c0,sim,bgc,time_series,forcing,MTM,PQ_inv), tol, parms);

% sol_fname = strcat(string(tName(sim.selection)),'_sol');
sol_fname = sprintf('%s_sol_ierr_%d', string(tName(sim.selection)), ierr);
fprintf('%s.m: Saving sol, it_hist, ierr" in %s"...\n', mfilename,sol_fname);
save(sol_fname, 'sol', 'it_hist','ierr','-v7.3','-nocompression');

myRestartFile = sprintf('%s/restart_%d_%s_x0.mat',             sim.outputRestartDir, round(sim.start_yr), strjoin(tName(sim.selection))      );
mySolFile     = sprintf('%s/restart_%d_%s_x0_sol_ierr_%d.mat', sim.outputRestartDir, round(sim.start_yr), strjoin(tName(sim.selection)), ierr);
fprintf('%s.m: Saving "%s" as the solution restart file "%s"...\n', mfilename,myRestartFile, mySolFile);
copyfile( myRestartFile, mySolFile);

num_r_iterations = 5;
% x = load('Data/restart_0_1_output_rIter/G_30.mat', 'x0').x0;
% load('Data_GP/saveRes_25.mat', 'x_hist');
% x = x_hist(:,25);
% clear x_hist

x = sol;
x_histx = x;
r_hist = zeros(size(x));
it_histx = zeros(num_r_iterations,1);
Npt = -1;


% num_r_iterations = 0;
fprintf('\n%s.m: Now we "relax" or forward integtration single variable solution a few iterations: %s\n',mfilename,datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
fprintf('%s.m: Starting Picard Integration for tracer selection = [%d]...\n',mfilename, sim.selection);
fprintf('%s.m: %d time steps, each time step is %1.1f (h), simulating %1.1f years\n', ...
    mfilename, sim.T/sim.dt, sim.dt/sim.const.sec_h, tot_t /sim.const.sec_y)

timer_loop = tic;
for itc = 1:num_r_iterations
    % FIXME: note "x" not "x0"

    [r,G, x1] = calc_G(x,c0,sim,bgc,time_series,forcing,MTM,PQ_inv);
    % disp('FIXME: hacking r = -G')
    % r = -G;
    % fnrm = norm(r);
    fprintf('%s.m: itc %d norm(G) = %g\n',mfilename, itc, norm(G))
    fprintf('%s.m: itc %d norm(r) = %g\n',mfilename, itc, norm(r))

    x_histx = [x_histx,x];
    r_hist = [r_hist,r];
    it_histx(itc,1) = norm(r);

    x = x1;

end


keyboard
sol = x;
for num_r_relax_iterations = 1:10
    fprintf("\n%s.m: starting relaxation year #%d\n", mfilename, num_r_relax_iterations)
    [sim, bgc, time_series] = phi(sim, bgc, time_series, forcing, MTM);
end
 
disp([mfilename,' finished...'])
elapsedTime_all_loc = toc(timer_loop);
disp(' ');
disp(['Runtime: ', num2str(elapsedTime_all_loc, '%1.0f'),' (s) or ', num2str(elapsedTime_all_loc/60, '%1.1f'), ' (m)'])
disp(['Runtime per location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps/sim.domain.num_wet_loc*1000, '%1.2f'), ' (ms) MARBL, advection, diffusion, mfactor()'])
disp(['Runtime all location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps, '%1.2f'),                    ' (s)  MARBL, advection, diffusion, mfactor()'])
disp(['Runtime all location per sim year : ', num2str(elapsedTime_all_loc/60/1440/tot_t*sim.const.sec_y, '%1.2f'), ' (d/y_sim)'])
disp(['Simulation speed: ', num2str(tot_t/elapsedTime_all_loc/sim.const.days_y, '%1.1f'), ' (sim y/d) aka (SYPD)'])

keyboard
% FIXME: need to save workspace?!
logDir = strcat(sim.outputRestartDir,'/Logs/');
if ~exist(logDir, 'dir')
       mkdir(logDir)
end
save_timer = tic; disp('Saving (possibly) large workspace file...'); save(strcat(logDir,'last_run.mat'),'-v7.3','-nocompression'); toc(save_timer);


disp(['Log file of one location for all time steps uses ',num2str(getMemSize(time_series)/1024/1024, '%1.1f'),' MB'])
disp(['Log file of one location uses ',num2str(getMemSize(time_series)/1024/sim.num_time_steps, '%1.1f'),' KB per time step'])
