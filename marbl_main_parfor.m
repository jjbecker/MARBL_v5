% function [singleColTracerError, singleColTracerSolution] = marbl_main_parfor(varargin) % tracer_loop, inputRestartFile, time_step_hr, logTracers, recalculate_PQ_inv, short_circuit

% "main" of cyclostationary transport version of MARBL,
% Try to solve MARBL using "Newton Like" methods from Kelley
%
%        https://ctk.math.ncsu.edu/newtony.html
%
%   matlab -nodisplay -nodesktop -nosplash -noFigureWindows -logfile batch.txt < marbl_nsoli.m &
%
% MARBL code being used is https://github.com/marbl-ecosys/MARBL

clear functions globals     % Need this to clear "persistent" variables in "G()", "time_step()" and "calculate_forcing()"
clearvars -except varargin;

timer_total = tic;
fprintf('%s.m: Start at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));

dbstop if error
format short g

addpath('MEX', genpath('utils'), genpath('plotting'));    % dname = sprintf('%s/../',myDataDir()); addpath(dname);

% Clean up for threads, more complex and much slower than expected.
% killAndCleanThreads();  % Leftover threads and dumps cause a lot of issues

diary off; diary off; diary on; diary off; diary on; diary on; % FIXME: diary behavior is such that if renamed it's still active diary!!

%%%%%%

% Setup big picture parts of a simulation and/or NK solution.

sim.runInParallel = 0;      % parfor can't use spmd inside, at least I can not make that work                       
sim.verbose_debug = 1;
sim.forwardIntegrationOnly    = 0;      % 1 -> no NK just fwd integration
sim.num_relax_iterations      = 0;      % 0 means no relax steps, just use NK x1_sol
sim.num_forward_years         = 0;      % if fwd only, num fwd, else this inum fwd after relax step

% setup file paths and selected tracers

% sim = setInputAndOutputFilePaths(sim, varargin)
sim = setInputAndOutputFilePaths(sim,[]);

% FIXME: hack in some stuff for debug

sim.tracer_loop = {'DOPr' 'DONr' 'Fe'};
sim.time_step_hr = 12;

sim.phi_years     = 1;      % NK always uses 1 year integration

%%%%%%
if 0 && sim.debug_disable_phi
    fprintf('\n\n\t%s.m: ********* phi() is short circuited skip MTM read  *********\n\n',mfilename)
    MTM = 1;
else
    load(sim.inputRestartFile,'tracer','state','MTM');
    % MTM        = load(sim.inputRestartFile,'MTM').MTM;
    bgc.tracer = tracer;    clear tracer;
    bgc.state  = state;     clear state;
end
% FIXME: do we?
% We need the initial tracers for the restart file[sim, 
%%%%% Need moles for calc_global_moles_and_means() for bgc2nsoli()
sim.grd     = load(sim.inputRestartFile,'sim').sim.grd;
sim.domain  = load(sim.inputRestartFile,'sim').sim.domain;
[sim, ~] = calculate_depth_map_and_volumes(sim);
sim = setPeek(sim);
sim = calc_global_moles_and_means(bgc, sim);

c0 = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
sz = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
% c  = reshape(c0, sz);
% x0_bgc  = replaceSelectedTracers(sim, c0, x0_sol, sim.selection);

%%%%%%% These are the output result from the "parfor slice variable", sized to accept all possible tracers
num_tr = numel(sim.tracer_loop);
singleColTracerError    = zeros([1, size(bgc.tracer,3)]);
singleColTracerSolution = zeros([size(sim.domain.iwet_FP,1), size(bgc.tracer,3)]);

%%%%%%
delete(gcp('nocreate'))
numCores = max(round(feature('numcores')/2), 6) % Never more than 6
numCores = min(num_tr, numCores)                % Never more than numTracer
if isunix && ismac
    numCores = min(2, numCores)     % hack for limited RAM on laptop
end
p = parpool(numCores);
if p.NumWorkers ~= numCores
    error("Could not allocate desired number of cores")
end
ticBytes(gcp);


% sim.recalculate_PQ_inv = 0;
sim
% keyboard
myFilename = 'marbl_main_parfor';     % mfilename does NOT work in a parfor
% parfor (par_idx = 1:num_tr,numCores)              % PARENTHESIS are CRUCIAL
for par_idx = 1:num_tr

    tmp_sim = sim;      % 'tmp_sim'  may be extraneous; parfor crazy!

    tracer_cell = tmp_sim.tracer_loop     (par_idx);
    tracer_col  = tmp_sim.tracer_loop_idx (par_idx);
    fprintf('%s.m: Starting tracer #%d (%s)\n', myFilename, tracer_col, string(tracer_cell))

%     [tmp_sim, tmp_bgc, tmp_ierr, tmp_x] = parfor_inner(tmp_sim, MTM, tracer_cell );
        [~, ~, ierr, x] = parfor_inner(tmp_sim, MTM, tracer_cell );
%         singleColTracerError   (:, tracer_col) = ierr;
%         singleColTracerSolution(:, tracer_col) = x;
end % of loop over tracers

tocBytes(gcp)
delete(gcp('nocreate'))

fprintf('...end of loop over tracers : '); toc(timer_total)
fprintf('Shutting down the parpool...\n')

bgc.tracer = nsoli2bgc(sim, bgc, x0_bgc);   % marbl format x0


% FIXME: need to save workspace?!
logDir = strcat(sim.outputRestartDir,'/Logs/');
if ~exist(logDir, 'dir')
    mkdir(logDir)
end
save_timer = tic; fprintf('Saving (possibly) large workspace file... '); save(strcat(logDir,'last_run.mat'),'-v7.3','-nocompression'); toc(save_timer);

fprintf('... end of %s.m ', mfilename);
elapsedTime_all_loops_all_tracers = toc(timer_total);
disp(' ');
fprintf('\n%s.m: Finished outer solution loops over %d tracers at %s\n', mfilename, numel(sim.tracer_loop),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
disp(['Runtime: ', num2str(elapsedTime_all_loops_all_tracers, '%1.0f'),' (s) or ', num2str(elapsedTime_all_loops_all_tracers/60, '%1.1f'), ' (m)'])
fprintf('%s.m: Finished at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));

% end % of function
