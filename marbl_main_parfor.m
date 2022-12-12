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
sim.phi_years     = 1;      % NK always uses 1 year integration


% FIXME: hack in some stuff for debug

sim.time_step_hr = 12;
% sim.recalculate_PQ_inv         = 0;
% sim.debug_disable_phi          = 1;
% sim.disable_ALL_Preconditioner = 1;


%%%%%%
sim.grd         = load(sim.inputRestartFile,'sim').sim.grd;
sim.domain      = load(sim.inputRestartFile,'sim').sim.domain;
[sim, ~]        = calculate_depth_map_and_volumes(sim);

sim  = setPeek(sim);

if 0 && sim.debug_disable_phi
    fprintf('\n\n\t%s.m: ********* phi() is short circuited skip MTM read  *********\n\n',mfilename)
    MTM = 1;
else
    fprintf('%s.m: Loading transports and tracers from restart file...\n', mfilename);
    load(sim.inputRestartFile,'tracer','state','MTM');
    % MTM        = load(sim.inputRestartFile,'MTM').MTM;
    bgc.tracer = tracer;    clear tracer;
    bgc.state  = state;     clear state;
end
num_tr_selected = numel(sim.tracer_loop);
sz_bgc          = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];

%%%%% 
% Need moles for calc_global_moles_and_means() for bgc2nsoli()
% We need initial tracers for restart file[sim, 
% FIXME: do we?
sim = calc_global_moles_and_means(bgc, sim);
c0  = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
c   = reshape(c0, sz_bgc);

%%%%%%
delete(gcp('nocreate'))
numCores = max(round(feature('numcores')/2), 6) % Never more than 6
numCores = min(num_tr_selected, numCores)                % Never more than numTracer
if isunix && ismac
    numCores = min(2, numCores)     % hack for limited RAM on laptop
end
p = parpool(numCores);
if p.NumWorkers ~= numCores
    error("Could not allocate desired number of cores")
end
ticBytes(gcp);


sim
keyboard
myFilename = 'marbl_main_parfor';     % mfilename does NOT work in a parfor
%%%%%%% 
% These are output result from "parfor slice variable", sized to accept all
% possible tracers.

singleColTracerError = zeros([1, size(bgc.tracer,3)]);
% singleColTracerSolution = 0 *c;
% FIXME: these temps DO work in loop below but end up inscrambled order...
tmp_xsol = zeros([size(sim.domain.iwet_FP,1), size(bgc.tracer,3)]);
tmp_ierr = zeros([1, size(bgc.tracer,3)]);

%%%

parforIdxRange = 1:num_tr_selected;

% par_idx is simply "order of execution" -NOT- tracer number
% so we nee3d tracerRange which is tracer number.
% 
% have to be SOOOOOO careful to avoid brodcasting huge struct. Using
% sim.tracer_loop will be a wasted copy of a huge sim struct.
tracerRange = sim.tracer_loop_idx;
tracer_cell = sim.tracer_loop;

% for par_idx = 1:3  % just for DEBUG
parfor (par_idx = parforIdxRange,numCores)  % PARENTHESIS are CRUCIAL

    % Also crucial: par_idx order randomly selected from range!

    % DEBUG have to be SOOOOOO careful to avoid brodcasting huge struct
    % tmp_sim = sim; %% NOT an issue since we have send whold sim anyway? 
    fprintf('%s.m: Starting tracer parfor idx %d which is for tracer %s (#%d)\n', myFilename, ...
        par_idx, string(tracer_cell(par_idx)), tracerRange (par_idx))

    % Works, but lots of code to catch two array, but easy to read...
    
    [~, ~, ierr, xsol] = parfor_inner(sim, MTM, string(tracer_cell(par_idx)));
    
    tmp_ierr (:, par_idx) = ierr;
    tmp_xsol (:, par_idx) = xsol;

    % Works, but hard to read...
    %
    % [~, ~, tmp_ierr(:,par_idx), tmp_xsol(:,par_idx)] = ...
    %     parfor_inner(sim, MTM, string(tracer_cell(par_idx)));
    % 
    %
    % FIXME: does NOT work, something about slicing. Most likely can only 
    % use parfor index, not one computed from it
    %
    % tracer_num  = tracerRange (par_idx);
    % singleColTracerError   (:, tracer_num) = ierr;
    % singleColTracerSolution(:, tracer_num) = x;
    
end % of loop over tracers

fprintf('...end of loop over tracers : '); toc(timer_total)
fprintf('Shutting down parpool...\n')

tocBytes(gcp)
delete(gcp('nocreate'))

%%%

% Now unscramble results captured in random order by parfor loop. No need
% for a for loop, use array index on slices.

tracerRange = sim.tracer_loop_idx (parforIdxRange);

singleColTracerError   (:, tracerRange) = tmp_ierr (:, parforIdxRange);
% singleColTracerSolution(:, tracerRange) = tmp_xsol (:, parforIdxRange);
c(:, tracerRange)                       = tmp_xsol (:, parforIdxRange); 

% Just replaced initial values in restart with single column solution. Now
% just need to create a restart file with this result.

keyboard
[~,oldName,~] = fileparts(sim.inputRestartFile)
newRestartFileName = sprintf('%s/%s_single_col_%s.mat', sim.outputRestartDir, oldName, strjoin(tracer_cell,'_'))
[sim, bgc] = saveRestartFiles(sim, bgc, c, newRestartFileName);

%%%%% whew. 

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
