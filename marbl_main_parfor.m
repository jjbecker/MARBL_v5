% MARBL code being used is https://github.com/marbl-ecosys/MARBL
% function [singleColTracerError, singleColTracerSolution] = marbl_main_parfor(varargin) % tracer_loop, inputRestartFile, time_step_hr, logTracers, recalculate_PQ_inv, short_circuit
% https://ctk.math.ncsu.edu/newtony.html
% matlab -nodisplay -nodesktop -nosplash -noFigureWindows -logfile batch.txt < marbl_nsoli.m &
%
% Clean up for threads, more complex and much slower than expected.
% killAndCleanThreads();    % Leftover threads and dumps cause a lot of issues
clear functions globals     % clear "persistent" vars
clearvars -except varargin;

format short g;
dbstop if error;
% dbstop if naninf;     % FIXME: intentional tracer NaNs triggers this
% dbstop if warning;    % FIXME: causes endless loop of obscure breaks

elapsedTimeTotal = tic;
fprintf('%s.m: Start at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));

addpath('MEX', genpath('utils'), genpath('plotting'));    % dname = sprintf('%s/../',myDataDir()); addpath(dname);
diary off; diary off; diary on; diary off; diary on; diary on; % FIXME: diary behavior is such that if renamed it's still active diary!!

%%%
% if outer loop is running use previous result...
clear newRestartFileName
for outerLoop_idx = 1:6

    clear calc_G calculate_forcing phi time_step_ann % clear "persistent" vars
    %%%
    % Setup big picture parts of a simulation.
    % sim = setInputAndOutputFilePaths([]);
    % sim = setInputAndOutputFilePaths(varargin)

    % tmpTracer_loop  = {'DOP' 'DOC'};
    tmpTracer_loop  = {'spFe' 'NH4'};
%     tmpTracer_loop  = tracer_names(0);
    tmpTime_step_hr = 3;

    tmpRecalculate_PQ_inv   = 1;    % default = 1
    tmpDebug_disable_phi    = 0;    % default = 0
    tmpLogTracer            = 0;    % default = 0

    % read an initial condition file.
    % assume for simplicity it is (probably) first pass...
    tmpInputFile = strcat(myDataDir(), 'restart_260_integrate_from_0.mat');
    % tmpInputFile = strcat(myDataDir(), 'restart_0_1_output/restart_260_integrate_from_0_DOP_DOC.mat');
    % tmpInputFile = strcat(myDataDir(), 'outerLoop_1.mat');
    %
    % outer loop #2 or greater?
    if exist('newRestartFileName','var')
        if tmpDebug_disable_phi
            fprintf('\n%s.m: outerLoop #%d initial condition WOULD HAVE BEEN %s but phi() is short circuited!! reuse initial condition  *********\n\n',mfilename, outerLoop_idx, newRestartFileName);
        else
            tmpInputFile = newRestartFileName;  % outer loop #2 or greater
        end
    end
    fprintf('%s.m: outerLoop #%d initial condition is %s\n', mfilename, outerLoop_idx, tmpInputFile);
    if ~exist(tmpInputFile, 'file')
        keyboard
    end

    % sim = setInputAndOutputFilePaths({tmpTracer_loop, tmpInputFile, tmpTime_step_hr});    % Order of args IMPORTANT !!!
    % sim = setInputAndOutputFilePaths({{'DOP' 'DOC'}, tmpInputFile, 12});                  % Order of args IMPORTANT !!!
    sim = setInputAndOutputFilePaths({tmpTracer_loop, tmpInputFile, tmpTime_step_hr, ...
        tmpRecalculate_PQ_inv, tmpDebug_disable_phi, tmpLogTracer });                       % Order of args IMPORTANT !!!
    % sim = setInputAndOutputFilePaths({{'DOP' 'DOC'}, tmpInputFile, 12, 0, 1, 0});         % Order of args IMPORTANT !!!

    clear tmpRecalculate_PQ_inv tmpDebug_disable_phi tmpLogTracer tmpTracer_loop tmpTime_step_hr tmpInputFile;

    %%%
    % Setup big picture parts of NK solution in marbl_solve():
    % relative tolerance; as fraction of f(x0)..
    sim.rtol     = 5e-1;    % marbl_solve(): stop if norm(drift,2) < 10% of G(x0)
    sim.maxfeval = 5;       % marbl_solve(): max number of function evaluation

    %%%
    % Never use more than numTracer
    numCores = min(round(feature('numcores')), numel(sim.tracer_loop));
    if (isunix && ismac)
        numCores  = min(2, numCores);   % limited RAM on laptop...
        numMatlab = numCores;
    else
        numCores = min(numCores, 10);  % be a good citizen on GP
        numMatlab= 2* numCores; % FIXME: does NOT work, only numCores
        numMatlab = numCores;
    end

    fprintf('%s.m: Start %d cores to run a total of %d Matlabs at %s\n', mfilename, numCores, numMatlab, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
    delete(gcp('nocreate'))
    p = parpool(numCores);
    if p.NumWorkers ~= numCores, error("Could not allocate desired number of cores"); end
    ticBytes(gcp);

    %%%
    % IMPORTANT     % have to be so careful to avoid brodcasting huge struct.
    % Using sim.tracer_loop will waste copy of a huge  struct.
    tracerRange = sim.tracer_loop_idx;
    tracer_cell = sim.tracer_loop;
    parforIdxRange = 1:numel(sim.tracer_loop);
    myFilename = mfilename;     % mfilename does NOT work in a parfor

    %%%
    [sim, bgc, MTM] = loadRestartFile(sim);
    sim

    % These are output result from parfor "slice variable", sized to accept all
    % possible tracers, in random order.
    tmp_xsol = zeros([size(sim.domain.iwet_JJ,1), size(bgc.tracer,3)]);
    tmp_ierr = zeros([1, size(bgc.tracer,3)]);
    tmp_fnrm = zeros([1, size(bgc.tracer,3)]);

    % eventually pass in result of first calculate of f(x0) which is same for
    % every tracer and save 3-4 hrs for every outer loop
    % f0 = calc_f0(sim, bgc, MTM, string(tracer_cell(par_idx)));

    fprintf('%s.m: "parfor" starting <=%d Matlab workers on %d cores\n', myFilename, numMatlab, numCores)
    % for par_idx = parforIdxRange  % DEBUG
    parfor (par_idx = parforIdxRange, numMatlab)  % PARENTHESIS are CRUCIAL

        % par_idx is (usually) randomly selected order from range!
        % par_idx is simply "order of execution" -NOT- tracer number
        % Need tracerRange which is tracer number "decoder"

        fprintf('%s.m: Starting tracer parfor idx %d which is for tracer %s (#%d)\n', myFilename, ...
            par_idx, string(tracer_cell(par_idx)), tracerRange (par_idx))

        % this code is easy to read and avoids crazy syntax for slice variable.

        [~, ~, ierr, sol, fnrm] = parfor_inner(sim, MTM, string(tracer_cell(par_idx)));

        tmp_ierr (:, par_idx) = ierr;
        tmp_xsol (:, par_idx) = sol;
        tmp_fnrm (:, par_idx) = fnrm;

% elapsedTimeSave = tic; 
% fprintf('%s.m: Finish integration of %s years: %s\n',mfilename,num2str(1+years_gone_by,2),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
%         fprintf('%s.m: Finished tracer parfor idx %d which is for tracer %s (#%d): %s\n', myFilename, ...
%             par_idx, string(tracer_cell(par_idx)), tracerRange (par_idx),datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')))
        fprintf('%s.m: Finished tracer parfor idx %d which is for tracer %s (#%d)\n', myFilename, ...
            par_idx, string(tracer_cell(par_idx)), tracerRange (par_idx))

    end % of loop over tracers

    fprintf('...end of loop over tracers : '); toc(elapsedTimeTotal)

    %%%
    % put results in correct order

    [newRestartFileName, singleColTracerError, singleColTracerNorm, c_sol]  = ...
        unscrambleParFor(outerLoop_idx, sim, bgc, tmp_ierr, tmp_xsol, tmp_fnrm, parforIdxRange, tracer_cell );

    %%% whew.
    fprintf('Shutting down parpool...\n')
    tocBytes(gcp)
    delete(gcp('nocreate'))

    logDir = strcat(sim.outputRestartDir,'/Logs/');
    if ~exist(logDir, 'dir')
        mkdir(logDir)
    end
    elapsedTimeSave = tic; fprintf('Saving (possibly) large workspace file... '); save(strcat(logDir,'last_run.mat'),'-v7.3','-nocompression'); toc(elapsedTimeSave);

    %%%
    fprintf('... end of %s.m ', mfilename);
    elapsedTimeAllLoopsAllTracers = toc(elapsedTimeTotal);
    disp(' ');
    fprintf('\n%s.m: Finished %d outer solution loops over %d tracers with as many as %d solution loops at %s\n', mfilename, outerLoop_idx, numel(sim.tracer_loop), sim.maxfeval, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
    disp([mfilename,'.m: Runtime: ', num2str(elapsedTimeAllLoopsAllTracers, '%1.0f'),' (s) or ', num2str(elapsedTimeAllLoopsAllTracers/60, '%1.1f'), ' (m)'])

end % of outerLoop_idx loop
disp(' ');
disp([mfilename,'.m: Runtime: ', num2str(elapsedTimeAllLoopsAllTracers, '%1.0f'),' (s) or ', num2str(elapsedTimeAllLoopsAllTracers/60, '%1.1f'), ' (m)'])
fprintf('%s.m: COMPLETELY Finished at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));


%%
function [sim, bgc, MTM] = loadRestartFile(sim)
if 0 && sim.debug_disable_phi
    fprintf('\n\n\t%s.m: ********* phi() is short circuited skip MTM read  *********\n\n',mfilename)
    MTM = 1;
else
    fprintf('%s.m: Loading transports and tracers from restart file...\n', mfilename);
    load(sim.inputRestartFile,'tracer','state','MTM');
    bgc.tracer = tracer;    clear tracer;
    bgc.state  = state;     clear state;
    % Need moles for calc_global_moles_and_means() in bgc2nsoli() or do we?
    sim = calc_global_moles_and_means(bgc, sim);
end
end % loadRestartFile

%%
function [newRestartFileName, singleColTracerError, singleColTracerNorm, c_sol]  = ...
    unscrambleParFor(outerLoop_idx, sim, bgc, tmp_ierr, tmp_xsol, tmp_fnrm, parforIdxRange, tracer_cell )

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

% [~,oldName,~] = fileparts(sim.inputRestartFile)
% newRestartFileName = sprintf('%s/%s_single_col_%s.mat', sim.outputRestartDir, oldName, strjoin(tracer_cell,'_'))
% newRestartFileName = sprintf('%s/%s_%s.mat', sim.outputRestartDir, oldName, strjoin(tracer_cell,'_'))
newRestartFileName = sprintf('%s/%s_%d.mat', sim.outputRestartDir, 'outerLoop', outerLoop_idx)

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

% %%
% function f0 = calc_f0(sim, bgc, MTM, tracer_str)
%
% [sim, bgc, bgc_struct, time_series, forcing] = init_sim(sim);
% % sim = calc_global_moles_and_means(bgc, sim);
%
% sz = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
% c0 = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
% c  = reshape(c0, sz);
%
% sim.selection = unique(sort(find( strcmp(tracer_names(0),tracer_str) )));
% x0 = c(:,sim.selection);    % initial condition for Nsoli()
% x0 = x0(:);                 % unitless
%
% PQ_inv = calc_PQ_inv(sim, bgc, time_series, forcing, MTM);
% f = @(x) calc_G(x, c0, sim, bgc, time_series, forcing, MTM, PQ_inv);
% f0=feval(f,x0);
%
% end % calc_f0
%
