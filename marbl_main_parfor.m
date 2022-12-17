% MARBL code being used is https://github.com/marbl-ecosys/MARBL
% function [singleColTracerError, singleColTracerSolution] = marbl_main_parfor(varargin) % tracer_loop, inputRestartFile, time_step_hr, logTracers, recalculate_PQ_inv, short_circuit
% https://ctk.math.ncsu.edu/newtony.html
% matlab -nodisplay -nodesktop -nosplash -noFigureWindows -logfile batch.txt < marbl_nsoli.m &

dbstop if error
format short g;
% fmt = format

% Clean up for threads, more complex and much slower than expected.
% killAndCleanThreads();    % Leftover threads and dumps cause a lot of issues
clear functions globals     % Need this to clear "persistent" variables in "G()", "time_step()" and "calculate_forcing()"
clearvars -except varargin;

timer_total = tic;
fprintf('%s.m: Start at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));

addpath('MEX', genpath('utils'), genpath('plotting'));    % dname = sprintf('%s/../',myDataDir()); addpath(dname);
diary off; diary off; diary on; diary off; diary on; diary on; % FIXME: diary behavior is such that if renamed it's still active diary!!

%%%
% Setup big picture parts of a simulation and/or NK solution.
% sim = setInputAndOutputFilePaths([]);
% sim = setInputAndOutputFilePaths(varargin)
% 
tmpTracer_loop  = {'DOP' 'DOC'};
tmpTracer_loop  = tracer_names(0);
tmpTime_step_hr = 3;
% tmpInputFile    = strcat(myDataDir(), 'Data_GP/restart_260_integrate_from_0.mat'); disp(tmpInputFile);
tmpInputFile    = strcat(myDataDir(), 'restart_260_integrate_from_0.mat'); disp(tmpInputFile);
% sim = setInputAndOutputFilePaths({tmpTracer_loop, tmpInputFile, tmpTime_step_hr});

% sim = setInputAndOutputFilePaths({{'DOP' 'DOC'}, tmpInputFile, 12});

tmpRecalculate_PQ_inv   = 1;
tmpDebug_disable_phi    = 0;
tmpLogTracer            = 0;
% Order of args IMPORTANT !!!
sim = setInputAndOutputFilePaths({tmpTracer_loop, tmpInputFile, tmpTime_step_hr, ...
    tmpRecalculate_PQ_inv, tmpDebug_disable_phi, tmpLogTracer });
% 
% sim = setInputAndOutputFilePaths({{'DOP' 'DOC'}, tmpInputFile, 12, 0, 1, 0});

%%%
% marbl_solve(): relative tolerance; as fraction of f(x0)..
sim.rtol     = 1e-2;    % marbl_solve(): stop if norm(drift,2) < 10% of G(x0)
sim.maxfeval = 2;       % marbl_solve(): max number of function evaluation

%%%
if 0 && sim.debug_disable_phi
    fprintf('\n\n\t%s.m: ********* phi() is short circuited skip MTM read  *********\n\n',mfilename)
    MTM = 1;
else
    fprintf('%s.m: Loading transports and tracers from restart file...\n', mfilename);
    load(sim.inputRestartFile,'tracer','state','MTM');
    bgc.tracer = tracer;    clear tracer;
    bgc.state  = state;     clear state;
end
% Need moles for calc_global_moles_and_means() in bgc2nsoli() or do we?
sim = calc_global_moles_and_means(bgc, sim);

%%%
% Never use more than numTracer or more than 1/2 nodes processors
numCores = min(round(feature('numcores')/2), numel(sim.tracer_loop)) 
if isunix && ismac, numCores = min(2, numCores) % limited RAM on laptop...
else numCores = min(numCores, 10)               % be a good citizen on GP
end

delete(gcp('nocreate'))
p = parpool(numCores);
if p.NumWorkers ~= numCores, error("Could not allocate desired number of cores"); end
ticBytes(gcp);

%%%
sim
% These are output result from parfor "slice variable", sized to accept all possible tracers, in random order.
tmp_xsol = zeros([size(sim.domain.iwet_JJ,1), size(bgc.tracer,3)]);
tmp_ierr = zeros([1, size(bgc.tracer,3)]);
tmp_fnrm = zeros([1, size(bgc.tracer,3)]);

% IMPORTANT     % have to be so careful to avoid brodcasting huge struct. 
                % Using sim.tracer_loop will waste copy of a huge  struct.
tracerRange = sim.tracer_loop_idx;
tracer_cell = sim.tracer_loop;
parforIdxRange = 1:numel(sim.tracer_loop);

myFilename = 'marbl_main_parfor';     % mfilename does NOT work in a parfor
% par_idx is simply "order of execution" -NOT- tracer number
% Need tracerRange which is tracer number "decoder"
% for par_idx = 1:1  % DEBUG
parfor (par_idx = parforIdxRange,numCores)  % PARENTHESIS are CRUCIAL

    % par_idx is (usually) randomly selected order from range!
    fprintf('%s.m: Starting tracer parfor idx %d which is for tracer %s (#%d)\n', myFilename, ...
        par_idx, string(tracer_cell(par_idx)), tracerRange (par_idx))

    % easy to read and avoids crazy syntax for slice variable.

    [~, ~, ierr, xsol, fnrm] = parfor_inner(sim, MTM, string(tracer_cell(par_idx)));

    tmp_ierr (:, par_idx) = ierr;
    tmp_xsol (:, par_idx) = xsol;
    tmp_fnrm (:, par_idx) = fnrm;

end % of loop over tracers

fprintf('...end of loop over tracers : '); toc(timer_total)

%%%
% put results in correct order

[singleColTracerError, singleColTracerNorm, c_sol]  = unscrambleParFor( sim, bgc, tmp_ierr, tmp_xsol, tmp_fnrm, parforIdxRange, tracer_cell );

%%% whew.
fprintf('Shutting down parpool...\n')
tocBytes(gcp)
delete(gcp('nocreate'))

logDir = strcat(sim.outputRestartDir,'/Logs/');
if ~exist(logDir, 'dir')
    mkdir(logDir)
end
save_timer = tic; fprintf('Saving (possibly) large workspace file... '); save(strcat(logDir,'last_run.mat'),'-v7.3','-nocompression'); toc(save_timer);

%%%
fprintf('... end of %s.m ', mfilename);
elapsedTime_all_loops_all_tracers = toc(timer_total);
disp(' ');
fprintf('\n%s.m: Finished outer solution loops over %d tracers at %s\n', mfilename, numel(sim.tracer_loop),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
disp(['Runtime: ', num2str(elapsedTime_all_loops_all_tracers, '%1.0f'),' (s) or ', num2str(elapsedTime_all_loops_all_tracers/60, '%1.1f'), ' (m)'])
fprintf('%s.m: Finished at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));

% end % of function


function [singleColTracerError, singleColTracerNorm, c_sol]  = ...
    unscrambleParFor( sim, bgc, tmp_ierr, tmp_xsol, tmp_fnrm, parforIdxRange, tracer_cell )

% Unscramble results captured in random order by parfor loop. No need
% for a for loop, use array index on slices. Then replace initial values in
% restart file with single column solution.


singleColTracerError = zeros([1, size(bgc.tracer,3)]);
singleColTracerNorm  = zeros([1, size(bgc.tracer,3)]);
% singleColTracerSolution = 0 *c0;

tracerRange = sim.tracer_loop_idx (parforIdxRange);
singleColTracerError (:, tracerRange) = tmp_ierr (:, parforIdxRange);
singleColTracerNorm  (:, tracerRange) = tmp_fnrm (:, parforIdxRange);
% tmp_c_sol          (:, tracerRange) = tmp_xsol (:, parforIdxRange);

ierrLimit = 1;

bad_par_idx = parforIdxRange( tmp_ierr >  ierrLimit )
badTracers  = sim.tracer_loop_idx ( bad_par_idx )

good_par_idx = parforIdxRange;
good_par_idx( bad_par_idx ) = []
goodTracers  = sim.tracer_loop_idx ( good_par_idx )

%%%
% Need moles for calc_global_moles_and_means() for bgc2nsoli()
% We need initial tracers for restart file[sim,
% FIXME: or do we?
c0_nsoli = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
sz_bgc = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
c_sol = reshape(c0_nsoli, sz_bgc);
c_sol (:, goodTracers) = tmp_xsol (:, good_par_idx);

[~,oldName,~] = fileparts(sim.inputRestartFile)
% newRestartFileName = sprintf('%s/%s_single_col_%s.mat', sim.outputRestartDir, oldName, strjoin(tracer_cell,'_'))
newRestartFileName = sprintf('%s/%s_%s.mat', sim.outputRestartDir, oldName, strjoin(tracer_cell,'_'))

tracer = nsoli2bgc(sim, bgc, c_sol);
[sim, bgc] = saveRestartFiles(sim, bgc, tracer, newRestartFileName);

tmp = reshape(c0_nsoli, sz_bgc);
for idx = parforIdxRange
    col = sim.tracer_loop_idx (idx);
    xnrm = norm(tmp(:,col)-c_sol(:,col));
    fprintf('%s.m: (%s)\tierr = %d fnrm(r) = %-#13.7g norm(sol-x0) = %-#10.7g\n', mfilename, string(tracer_cell(idx)), ...
        singleColTracerError(col), singleColTracerNorm(col), xnrm)
end

end % unscrambleParFor function


% 
% singleColTracerError = zeros([1, size(bgc.tracer,3)]); 
% singleColTracerNorm  = zeros([1, size(bgc.tracer,3)]);
% % singleColTracerSolution = 0 *c0;
% 
% tracerRange = sim.tracer_loop_idx (parforIdxRange);
% singleColTracerError (:, tracerRange) = tmp_ierr (:, parforIdxRange);
% singleColTracerNorm  (:, tracerRange) = tmp_fnrm (:, parforIdxRange);
% % tmp_c_sol          (:, tracerRange) = tmp_xsol (:, parforIdxRange);
% 
% tracerRange = sim.tracer_loop_idx (parforIdxRange);
% 
% ierrLimit = 1;
% 
% bad_par_idx = parforIdxRange( tmp_ierr >  ierrLimit )
% badTracers  = sim.tracer_loop_idx ( bad_par_idx )
% 
% good_par_idx = parforIdxRange;
% good_par_idx( bad_par_idx ) = []
% goodTracers  = sim.tracer_loop_idx ( good_par_idx )
% 
% c_sol = reshape(c0_nsoli, sz_bgc);
% c_sol (:, goodTracers) = tmp_xsol (:, good_par_idx);
% 
% tmp = reshape(c0_nsoli, sz_bgc);
% for idx = parforIdxRange
%     col = sim.tracer_loop_idx (idx);
%     xnrm = norm(tmp(:,col)-c_sol(:,col));
%     fprintf('%s.m: (%s)\tierr = %d fnrm(r) = %-#13.7g norm(sol-x0) = %-#10.7g\n', myFilename, string(tracer_cell(idx)), ...
%         singleColTracerError(col), singleColTracerNorm(col), xnrm)
% end
% still need to create a restart file with this result.
% keyboard
% 
% [~,oldName,~] = fileparts(sim.inputRestartFile)
% % newRestartFileName = sprintf('%s/%s_single_col_%s.mat', sim.outputRestartDir, oldName, strjoin(tracer_cell,'_'))
% newRestartFileName = sprintf('%s/%s_%s.mat', sim.outputRestartDir, oldName, strjoin(tracer_cell,'_'))
% 
% tracer = nsoli2bgc(sim, bgc, c_sol);
% [sim, bgc] = saveRestartFiles(sim, bgc, tracer, newRestartFileName);
% 
