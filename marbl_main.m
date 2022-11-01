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

clear all; % Need this to clear "persistent" variables in "G()", "time_step()" and "calculate_forcing()"
timer_total = tic;

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

%%%%%%

% Setup big picture parts of a simulation and/or NK solution.

verbose_debug = 1;

% Input restart file

% start_yr = 4101;  inputRestartFileStem = 'Data/InputFromAnn/restart4101.mat';
% start_yr =   0;   inputRestartFileStem = 'Data/passive_restart_init.mat'; % from netCDF 5/25/22
% start_yr = 260;   inputRestartFileStem = 'Data_GP/restart_260_integrate_from_0.mat';
start_yr = 1323;  inputRestartFileStem = 'restart_0_1_output/restart_1323_DOP_sol_x1.mat';
% clear start_yr

inputRestartFile = strcat(myDataDir(), inputRestartFileStem);
%     [outputArg1,outputArg2] = setInputAndOutputFilePaths(inputArg1,inputArg2);
fprintf('%s.m: Loop over tracers starts with: %s\n', mfilename, inputRestartFile);

% always need a selected tracer! For plot time series, or solve
% Most of these work very well in single tracer solution...
% FIXME: rather than get accurate solution, loop over tracers and try to reduce G
% to 1% of starting value, while looping over all tracers, and then
% repeat until to get final very accurate result where G is sqrt(eps)

% tracer_loop = {'DOPr' 'DONr' 'DOCr' 'O2' 'DON' 'DOC' 'DOP' 'diatSi' 'spCaCO3' 'diazFe' };
% tracer_loop = {'DOC' 'DOP' 'spCaCO3' 'diatSi' 'diazFe' };
tracer_loop = {'spCaCO3' 'diatSi' 'diazFe' };
tracer_loop = {'O2' };

tName = tracer_names(0);    % no CISO tracers
for tracer_str = tracer_loop
    % Need this to clear "persistent" variables in "G()"
    clear calc_G

    % matches(tName,tracer_str)
    selection = [ find( strcmp(tName,tracer_str) ) ];

    forwardIntegrationOnly = 0; % 1 -> no NK just fwd integration
    recalculate_PQ_inv     = 1; % recalculate J, PQ,PQ_inv or load file

    % remember brsola() "sol" is x0 value. x1 value is NOT last col of
    % x_hist; it is sol !!!
    % First relax iteration of x0 gives same x1 as sol run.
    % To be useful num_relax_iterations >= 2 if using x0_sol, but OK for x1_sol

    num_relax_iterations      = 0;      % 0 means no relax steps, just use NK x1_sol
    num_forward_years         = 0;      % if fwd only, num fwd, else this inum fwd after relax step
    yearsBetweenRestartFiles  = 10;
    logTracers                = 1;
    captureAllSelectedTracers = 0;

    % FIXME: Someday, when we know what inputs need to be, put all this a file
    time_step_hr              = 3;
    phi_years                 = 1;      % NK always using 1 year integration
    debug_PQ_inv              = 0;
    debug_disable_phi         = 0;

    %%%%
    % DEBUG stuff
% logTracers         = 0;
time_step_hr       = 12; % FAST debug
% debug_PQ_inv       = 1
% debug_disable_phi  = 1
recalculate_PQ_inv = 0

    %%%%%%
    marbl_file = 'Data/marbl_in'; % MARBL chemistry and other constants.
    %%%%%%

    fprintf('%s.m: Solving for tracer: %s\n', mfilename, string(tracer_str));
    fprintf('%s.m: Reading (Matlab) OFFLINE sim restart file with tracers and transports: %s\n', mfilename, inputRestartFile);
    % load() does NOT need file extension, but copy() does. sigh

    if ~isfile(inputRestartFile)
        error("missing file or typo in name of inputRestartFile")
    end
    load(inputRestartFile,'sim','MTM');

    % We just over wrote sim struct, so now we can save stuff in it again. sigh

    sim.verbose_debug           = verbose_debug;
    sim.forwardIntegrationOnly  = forwardIntegrationOnly ;
    sim.inputRestartFile        = inputRestartFile;
    sim.start_yr                = start_yr;
    sim.selection               = selection;
    sim.captureAllSelectedTracers=captureAllSelectedTracers;
    sim.logTracers              = logTracers;
    sim.logDiags                = and (0, sim.logTracers) ; % Usually no diags..
    sim.checkNeg = 0;


    sim.debug_PQ_inv            = debug_PQ_inv;
    sim.debug_disable_phi       = debug_disable_phi;


    % FIXME: lots of old and or leftover broken code floating around, clear tmp vars...
    clear inputRestartFile
    clear selection captureAllSelectedTracers logTracers

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

    sim.tot_t = sim.dt*sim.num_time_steps;  % used only for debug output

    % DEBUG: run on a few locations, rather global, MUCH faster debug something
    % sim.domain.num_wet_loc = 1; % comment out too run entire world

    % scale units to unity, using values of global mean of initial
    % need a column vector

    sim = calc_global_moles_and_means(bgc, sim);

    clear forwardIntegrationOnly marbl_file

    %%%%%%

    % =============== This is NK solver code ================

    toc(timer_total)

    fprintf('\n%s.m: Start Newton (Broyden Method) solver: %s\n',mfilename,datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
    timer_PQ_init_solve_relax_fwd = tic;


    % NK always using 1 year integration

    if phi_years ~= 1
        disp('ERROR: NK requires phi_years = 1')
        keyboard
    end

    %%%%
    % ALWAYS punt "ALT" methods that do NOT depend on any other tracers...

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
        if(recalculate_PQ_inv)
            PQ_inv = calc_PQ_inv(sim, bgc, time_series, forcing, MTM);
        else
            fprintf('\n%s.m: Loading ~30 GB(!) mfactored preconditioner PQ_inv from %s solution...\n', mfilename, strcat(string(tName(sim.selection))))
            tStart = tic;
            if sim.debug_PQ_inv
                PQ_inv = 1
            else
                load (strcat(myDataDir(),'sol/',strjoin(tName(sim.selection)),'_QJ'), 'PQ_inv')
            end
            fprintf('%s.m: %1.0f (s) to init sim and load PQinv \n',mfilename, toc(tStart));
        end % calculate or load PQ_inv



        [ierr, myRestartFile_x0, x0_sol, c0, sim, bgc, time_series, forcing, MTM, PQ_inv] = ...
            marbl_solve(x0, c0, sim, bgc, time_series, forcing, MTM, PQ_inv);


        % FIXME: use x0 or x1 of solution?
        x = x0_sol;

        if num_relax_iterations > 0
            [x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv, myRestartFile_relaxed] = ...
                marbl_relax(num_relax_iterations, x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv);
        end % relax step

    end % Solve for selected tracer

    % Next! allow ALL tracers, not just selection, to "relax' to solution.
    %    do pure forward integration for a while...

    years_gone_by = 0;
    for fwd_itc = 1:num_forward_years

        fprintf("\n%s.m: starting forward integrate year #%d of %d\n", mfilename, fwd_itc, num_forward_years)

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
    disp([mfilename,' finished ', strjoin(tName(sim.selection))])
    disp(['Runtime: ', num2str(elapsedTime_all_loc, '%1.0f'),' (s) or ', num2str(elapsedTime_all_loc/60, '%1.1f'), ' (m)'])
    %     disp(['Runtime per location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps/sim.domain.num_wet_loc*1000, '%1.2f'), ' (ms) MARBL, advection, diffusion, mfactor()'])
    %     disp(['Runtime all location per iteration: ', num2str(elapsedTime_all_loc/sim.num_time_steps, '%1.2f'),                    ' (s)  MARBL, advection, diffusion, mfactor()'])
    %     disp(['Runtime all location per sim year : ', num2str(elapsedTime_all_loc/60/1440/sim.tot_t*sim.const.sec_y, '%1.2f'), ' (d/y_sim)'])
    disp(['Simulation speed: ', num2str(sim.tot_t/elapsedTime_all_loc/sim.const.days_y, '%1.1f'), ' (sim y/d) aka (SYPD)'])
    disp(' ');

    inputRestartFile = myRestartFile_x0;

end % of loop over tracers
fprintf('...end of loop over tracers : '); toc(timer_total)


% FIXME: need to save workspace?!
logDir = strcat(sim.outputRestartDir,'/Logs/');
if ~exist(logDir, 'dir')
    mkdir(logDir)
end
save_timer = tic; fprintf('Saving (possibly) large workspace file... '); save(strcat(logDir,'last_run.mat'),'-v7.3','-nocompression'); toc(save_timer);

fprintf('... end of %s.m ', mfilename); toc(timer_total)
