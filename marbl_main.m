function [ierr, x0_sol] = marbl_main(varargin) % tracer_loop, inputRestartFile, time_step_hr, logTracers, recalculate_PQ_inv, short_circuit
% function varargout = redplot(varargin)
%     [varargout{1:nargout}] = plot(varargin{:},'Color',[1,0,0]);
% end


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

clear functions globals
clearvars -except varargin;

% Need this to clear "persistent" variables in "G()", "time_step()" and "calculate_forcing()"
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

forwardIntegrationOnly    = 0; % 1 -> no NK just fwd integration
num_relax_iterations      = 0;      % 0 means no relax steps, just use NK x1_sol
num_forward_years         = 0;      % if fwd only, num fwd, else this inum fwd after relax step

sim.verbose_debug = 1;
sim = setInputAndOutputFilePaths(sim, varargin)

tName = tracer_names(0);    % no CISO tracers
if ~all(matches(sim.tracer_loop,tName))
    %         error('%s.m: tracer %s not in {%s}', mfilename, string(tracer_str), strjoin(tName));
    error('\n%s.m: tracer list "%s" contains invalid tracer name', mfilename, strjoin(string(sim.tracer_loop)));
end


for tracer_str = sim.tracer_loop

    fprintf('%s.m: Solving for tracer: %s\n', mfilename, string(tracer_str));
    % ALWAYS punt "ALT" methods that do NOT depend on any other tracers...
    % remove duplicates and make sure all the choices are valid...
    sim.selection = [ find( strcmp(tName,tracer_str) ) ];
    sim.selection(ismember(sim.selection, [9,11]))=[];
    sim.selection = unique(sort(sim.selection));
    cstr = tName(sim.selection)';
    fprintf('%s.m: Selected tracers: %s \n', mfilename, cstr{:});


    % Need this to clear "persistent" variables in "G()" and "phi()"
    clear calc_G phi

    fprintf('%s.m: Reading (Matlab) OFFLINE sim restart file with tracers and transports: %s\n', mfilename, sim.inputRestartFile);
    if ~isfile(sim.inputRestartFile)
        error("missing file or typo in name of inputRestartFile")
    end

    MTM = load(sim.inputRestartFile,'MTM').MTM;
    sim.grd     = load(sim.inputRestartFile,'sim').sim.grd;
    sim.domain  = load(sim.inputRestartFile,'sim').sim.domain;

    sim = setPeek(sim);

    sim.phi_years = 1;      % NK always uses 1 year integration

    % remember brsola() "sol" is x0 value. x1 value is NOT last col of
    % x_hist; it is sol !!!
    % First relax iteration of x0 gives same x1 as sol run.
    % To be useful num_relax_iterations >= 2 if using x0_sol, but OK for x1_sol
    %%%%%% End of "inputs"

    % Stuff below is not a simulation input. It is code to setup grids, etc a
    % forward integration, or NK(), or...
    %

    %%%%%%
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


    if( forwardIntegrationOnly )
        fprintf('%s.m: forward integration ONLY\n',mfilename);
    else
        % Solve for selected tracer
        if(sim.recalculate_PQ_inv)
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

    sim.inputRestartFile = myRestartFile_x0;

end % of loop over tracers
fprintf('...end of loop over tracers : '); toc(timer_total)


% FIXME: need to save workspace?!
logDir = strcat(sim.outputRestartDir,'/Logs/');
if ~exist(logDir, 'dir')
    mkdir(logDir)
end
save_timer = tic; fprintf('Saving (possibly) large workspace file... '); save(strcat(logDir,'last_run.mat'),'-v7.3','-nocompression'); toc(save_timer);

fprintf('... end of %s.m ', mfilename); toc(timer_total)

end

    %     % FIXME: Someday, when we know what inputs need to be, put all this a file
    %     time_step_hr              = 3;
    %     debug_PQ_inv              = 0;
    %     debug_disable_phi         = 0;
    %
    %     %%%%
    %     % DEBUG stuff
    % % time_step_hr       = 12; % FAST debug
    % logTracers         = 0;
    % recalculate_PQ_inv = 0
    % debug_PQ_inv       = 1
    % debug_disable_phi  = 1

    %     sim.verbose_debug           = verbose_debug;
    %     forwardIntegrationOnly  = forwardIntegrationOnly ;
    %     sim.inputRestartFile        = inputRestartFile;
    %     sim.start_yr                = start_yr;
    %     sim.selection               = selection;
    %     sim.captureAllSelectedTracers=captureAllSelectedTracers;
    %     sim.logTracers              = logTracers;
    %     sim.checkNeg = 0;
    %
    %
    %     sim.debug_PQ_inv            = debug_PQ_inv;
    %     sim.debug_disable_phi       = debug_disable_phi;
    %
    %
    %     % FIXME: lots of old and or leftover broken code floating around, clear tmp vars...
    %     clear inputRestartFile
    %     clear selection captureAllSelectedTracers logTracers

