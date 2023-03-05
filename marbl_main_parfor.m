% MARBL code being used is https://github.com/marbl-ecosys/MARBL
% function [tracerError,x,y,z tracerSolution] = marbl_main_parfor(varargin) % tracer_loop, inputRestartFile, time_step_hr, logTracers, recalculate_PQ_inv, short_circuit
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
clear newRestartFileName tmpTracer_loop

% Solve just the "inorganic" tracers, let the biology settle in the "num_single_tracer_relax_iters" loop
% tmpTracer_loop  = tracer_names(0);
% ignore biology in solution, solves implicity with time steps after
% solving PO4 etc...
tName = tracer_names(0);
% tmpTracer_loop  = tName([1:32]);
tmpTracer_loop  = tName([1:17]);
% tmpTracer_loop  = tName([18:32]);

% tmpTracer_loop = {'PO4' 'NO3' 'SiO3' 'NH4' 'Fe' 'Lig' 'O2' 'DIC' 'DIC_ALT_CO2' 'ALK' 'ALK_ALT_CO2' 'DOC' 'DON' 'DOP' 'DOPr' 'DONr' 'DOCr'};
% DIC ALK SiO3 already stable, just wastes time to update. 
% PO4 NO3 SiO3 might be stable, just wastes time to update. 
% Fe not stable, and does not solve single tracer, just wastes time to update. 
% PO4 NO3 SiO3 might be stable, but do not solve single tracer, just wastes time to update. 
% spCaCO3 clearly diverges if not single tracer solved
% diazFe might diverge if not single tracer solved
% tmpTracer_loop  = {'Fe' 'DIC'}
% tmpTracer_loop  = {'O2'};

% ignore "DIC_ALT" and "ALK_ALT"
tmpTracer_loop(ismember(tmpTracer_loop,{'DIC_ALT_CO2' 'ALK_ALT_CO2'}) >0) = [];

% Fe not stable, and does not solve single tracer, just wastes time to update.
% "" NH4 DIC ALK DOP spChl diatChl diazChl
% tmpTracer_loop(ismember(tmpTracer_loop,{'NH4' 'Fe' 'DIC' 'ALK' 'DOP' 'spChl' 'diatChl' 'diazChl' 'spC' 'diatC' }) >0) = [];
tmpTracer_loop(ismember(tmpTracer_loop,{  'NH4' 'Fe' 'DIC' 'ALK' 'DOP' 'spChl' 'diatChl' 'diazChl' }) >0) = [];

% Shuffle tracers: TRY to avoid blocking by slower tracers in parfor loop.
tmpTracer_loop = tmpTracer_loop ( randperm ( length ( tmpTracer_loop )))

numOuterLoops = 10;
% numOuterLoops = 2;
for outerLoop_idx = 1:numOuterLoops

    clear calc_G calculate_forcing phi time_step_ann % clear "persistent" vars
    %%%
    % Setup big picture parts of a simulation.
    % sim = setInputAndOutputFilePaths([]);
    % sim = setInputAndOutputFilePaths(varargin)

    tmpTime_step_hr = 3;
% tmpTime_step_hr = 12;

    tmpRecalculate_PQ_inv   = 1;    % default = 1
    tmpDebug_disable_phi    = 0;    % default = 0
    tmpLogTracer            = 1;    % default = 1
% tmpRecalculate_PQ_inv   = 0;    % default = 1
% tmpDebug_disable_phi    = 1;    % default = 0
% tmpLogTracer            = 1;    % default = 0

    % read an initial condition file.
    % FIXME: Matlab can NOT use chmod a+w "locked" attribute set in the Finder. Have to make a writtable copy of inputRestartFile 
    % assume for simplicity it is (probably) first pass...
    % tmpInputFile = strcat(myDataDir(), 'restart_260_integrate_from_0.mat');
    % tmpInputFile = strcat(myDataDir(), 'restart_0_1_output/restart_260_integrate_from_0_DOP_DOC.mat');
    % tmpInputFile = strcat(myDataDir(), 'outerLoop_0.mat');
    % tmpInputFile = strcat(myDataDir(), 'tmp/outerLoop_1_tmp.mat');    % restart from outerLoop_1_tmp
    tmpInputFile = strcat(myDataDir(), 'passive_restart_init.mat');
    
    % outer loop #2 or greater?
    if exist('newRestartFileName','var')
        %         if tmpDebug_disable_phi
        %             fprintf('\n%s.m: outerLoop #%d initial condition WOULD HAVE BEEN %s but phi() is short circuited!! reuse initial condition  *********\n\n',mfilename, outerLoop_idx, newRestartFileName);
        %         else
        tmpInputFile = newRestartFileName;  % outer loop #2 or greater
        %         end
    end
    fprintf('%s.m: outerLoop #%d initial condition is %s\n', mfilename, outerLoop_idx, tmpInputFile);
    if ~exist(tmpInputFile, 'file')
        keyboard
    end

    % sim = setInputAndOutputFilePaths({tmpTracer_loop, tmpInputFile, tmpTime_step_hr});    % Order of args IMPORTANT !!!
    % sim = setInputAndOutputFilePaths({{'DOP' 'DOC'}, tmpInputFile, 12});                  % Order of args IMPORTANT !!!
    % sim = setInputAndOutputFilePaths({{'DOP' 'DOC'}, tmpInputFile, 12, 0, 1, 0});         % Order of args IMPORTANT !!!
    sim = setInputAndOutputFilePaths({tmpTracer_loop, tmpInputFile, tmpTime_step_hr, ...
        tmpRecalculate_PQ_inv, tmpDebug_disable_phi, tmpLogTracer });                       % Order of args IMPORTANT !!!

    clear tmpRecalculate_PQ_inv tmpDebug_disable_phi tmpLogTracer tmpTime_step_hr tmpInputFile;

    %%%
    % Setup big picture parts of NK solution in marbl_solve():
    % relative tolerance; as fraction of f(x0)..
sim.rtol     = 5e-1;        % marbl_solve(): stop if norm(drift,2) < 10% of G(x0)
% sim.maxfeval = 3;           % marbl_solve(): max number of function evaluation
sim.maxfeval = 10;           % marbl_solve(): max number of function evaluation
sim.num_forward_iters = 3;  % years of all tracer relax; aka num of bgc = phi(bgc) loops after marbl_solve.
% sim.maxfeval = 1;           % DEBUG ONLY
% sim.num_forward_iters = 1;  % DEBUG ONLY

    %%%
    % Set some reasonable limit...
    % Never use more than half of cores available for workers on GP.
    % total cores = 1 client + as many workers as we can get away with

    if ismac
        maxCores = ceil(feature('numcores'))-2; % save one for user
        maxCores = 2; % laptop only has 32 GB of RAM, jobs need 30/tracer...
    else
        maxCores = ceil(feature('numcores')/2); % Green Planet cluster
    end
%     maxCores = maxCores -1;                     % one for client
    if numel(sim.tracer_loop) <= maxCores       % small job just run it
        numCores = numel(sim.tracer_loop);   
    elseif ( ceil(numel(sim.tracer_loop)/2) <= maxCores ) % ad hoc; 2 batchs
        numCores = ceil(numel(sim.tracer_loop)/2);    
    else
        numCores = maxCores;
    end
    fprintf('numCores = %d\tnumTracers = %d\t workerCores = %d\n',  feature('numcores'), numel(sim.tracer_loop), numCores )
    
    fprintf('%s.m: "parfor" using %d workers +1 client = %d of %d cores on %d tracers...\n', mfilename, numCores, numCores+1, feature('numcores'), numel(sim.tracer_loop))
    parforIdxRange = 1:numel(sim.tracer_loop);
    delete(gcp('nocreate'))
    p = parpool(numCores);
    % parfor might start, but NOT get requested parpool
    if p.NumWorkers ~= numCores, error("Could not allocate desired number of cores"); end
%     ticBytes(gcp);

    %%%
    myFilename = mfilename;     % mfilename does NOT work in a parfor
    % IMPORTANT     % have to be so careful to avoid brodcasting huge struct.
    % Using sim.tracer_loop will waste copy of a huge  struct.
    tracerRange = sim.tracer_loop_idx;
    tracer_cell = sim.tracer_loop;

    %%%
    [sim, bgc, MTM] = loadRestartFile(sim);
    min(bgc.tracer(:))                           % force nonnegative
    bgc.tracer = max(-sqrt(eps), bgc.tracer);    % force nonnegative

    sim
    sim.tracer_loop
    %%%
    % These are output result from parfor "slice variable", sized to accept
    % all possible tracers, in random order
    tmp_xsol = zeros([size(sim.domain.iwet_JJ,1), size(bgc.tracer,3)]);
    tmp_ierr = zeros([1, size(bgc.tracer,3)]);
    tmp_fnrm = zeros([1, size(bgc.tracer,3)]);

    % for par_idx = parforIdxRange  % DEBUG
    % for par_idx = 1:0                             % restart from outerLoop_1_tmp
    % parfor (par_idx = parforIdxRange, numMatlab)  % PARENTHESIS are CRUCIAL
    parfor (par_idx = parforIdxRange)             % PARENTHESIS are CRUCIAL

        % par_idx is (usually) randomly selected order from range!
        % par_idx is simply "order of execution" -NOT- tracer number
        % Need tracerRange which is tracer number "decoder"

        fprintf('%s.m: Starting tracer parfor idx %d which is for tracer %s (#%d)\n', myFilename, ...
            par_idx, string(tracer_cell(par_idx)), tracerRange (par_idx))

        % this code is easy to read and avoids crazy syntax for slice variable.
        % do NOT need tracer values in sim, bgc, or time_series, but do
        % need stucts definition and so on later for forward/relax run

        [~, ~, ~, ~, ierr, sol , fnrm] = parfor_inner(sim, MTM, string(tracer_cell(par_idx)));

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
%     tocBytes(gcp)
    fprintf('Shutting down parpool...\n')
    delete(gcp('nocreate'))     % close parpool parfor created...

    %%%
    % put results in correct order
    % result is saved in a file...

    % newRestartFileName = sprintf('%s/%s_%d_tmp.mat', sim.outputRestartDir, 'outerLoop', outerLoop_idx);
    newRestartFileName = sprintf('%s/%s_%d.mat', sim.outputRestartDir, 'outerLoop', outerLoop_idx);
    if sim.num_forward_iters >0
        newRestartFileName = sprintf('%s/%s_%d_tmp.mat', sim.outputRestartDir, 'outerLoop', outerLoop_idx);
    end
    if sim.debug_disable_phi
        fprintf('\n\n\t%s.m: ********* phi() is short circuited; just copy newRestartFileName  *********\n\n',mfilename);
% FIXME: Matlab can NOT use chmod a+w "locked" attribute set in the Finder. Have to make a writtable copy of inputRestartFile
        copyfile( sim.inputRestartFile, newRestartFileName); % FIXME: this can NOT over ride a "locked file" set in the Finder. 
        fileattrib( sim.inputRestartFile,'+w','a');
    end % if

    [tracerError, tracerNorm, c_sol]  = unscrambleAndSaveParFor(newRestartFileName, sim, bgc, tmp_ierr, tmp_xsol, tmp_fnrm, parforIdxRange, tracer_cell );


    %%%
    % POSSIBLY allow ALL tracers, not just selection, to "relax' to
    % solution by doing pure forward integration for a while.

    if sim.num_forward_iters >0

        tName = tracer_names(0);    % no CISO tracers
        %sim.selection is only selects debug plots for desired tracer
        sim.selection = unique(sort(find( strcmp(tName,string(tracer_cell(parforIdxRange(1)))) )));

        % all we need to do is suck combined tracer from combined file and
        % use forcing and time_seried from parfor_inner, but can NOT return
        % sim and bgc so go thru all these gyrations...
        sim.phi_years = sim.num_forward_iters;
        % sim.inputRestartFile = newRestartFileName;          % restart from outerLoop_1_tmp by NOT using the combo which is fake.
        [sim, bgc, ~, time_series, forcing] = init_sim(sim);
        % init_sim reads the initial condition file so force nonnegative -and- restart from outerLoop_1_tmp by
        min(bgc.tracer(:))                           % force nonnegative -and- restart from outerLoop_1_tmp by
        bgc.tracer = max(-sqrt(eps), bgc.tracer);    % force nonnegative -and- restart from outerLoop_1_tmp by

        % now we can run forward a few years to couple the tracers we did
        % NOT solve with nsoli()
        
        [sim, bgc, time_series] = phi(sim, bgc, time_series, forcing, MTM);
        
        sim.start_yr  = sim.start_yr+1;

        % be sure to shutdown MEX
        mex_marbl_driver('shutdown');

        % save forward iteration result...
        newRestartFileName = sprintf('%s/%s_%d.mat', sim.outputRestartDir, 'outerLoop', outerLoop_idx);
        saveRestartFiles(sim, bgc.tracer, newRestartFileName);

        % POSSIBLY need to reset time span of a phi() first
        sim.phi_years = 1;

    end % if

    %%% whew.
    fprintf('Shutting down parpool...\n')
%     tocBytes(gcp)
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
function [tracerError, tracerNorm, c_sol]  = unscrambleAndSaveParFor(newRestartFileName, sim, bgc, tmp_ierr, tmp_xsol, tmp_fnrm, parforIdxRange, tracer_cell )

% Unscramble results captured in random order by parfor loop. No need
% for a for loop, use array index on slices. Then replace initial values in
% restart file with single column solution.

tracerError = zeros([1, size(bgc.tracer,3)]);
tracerNorm  = zeros([1, size(bgc.tracer,3)]);
% singleColTracerSolution = 0 *c0;

tracerRange = sim.tracer_loop_idx (parforIdxRange);
tracerError (:, tracerRange) = tmp_ierr (:, parforIdxRange);
tracerNorm  (:, tracerRange) = tmp_fnrm (:, parforIdxRange);
% tmp_c_sol          (:, tracerRange) = tmp_xsol (:, parforIdxRange);

ierrLimit = 1;

bad_par_idx = parforIdxRange( tmp_ierr >  ierrLimit )
badTracers  = sim.tracer_loop_idx ( bad_par_idx )

good_par_idx = parforIdxRange;
good_par_idx( bad_par_idx ) = []
goodTracers  = sim.tracer_loop_idx ( good_par_idx )

%%%
% Need moles for calc_global_moles_and_means() for bgc2nsoli()
% Need initial tracers for restart file

c0_nsoli = bgc2nsoli(sim, bgc.tracer);    % nsoli format; unitless; aka scaled FP
sz_bgc = [ numel(sim.domain.iwet_JJ) , size(bgc.tracer,3) ];
c_sol = reshape(c0_nsoli, sz_bgc);
c_sol (:, goodTracers) = tmp_xsol (:, good_par_idx);

tracer = nsoli2bgc(sim, bgc, c_sol);
% result is saved in a file, bgc and sim returned are bogus...
saveRestartFiles(sim, tracer, newRestartFileName);

% debug output
tmp = reshape(c0_nsoli, sz_bgc);
for idx = parforIdxRange
    col = sim.tracer_loop_idx (idx);
    xnrm = norm(tmp(:,col)-c_sol(:,col));
    fprintf('%s.m: (%s)\tierr = %d fnrm(r) = %-#13.7g norm(sol-x0) = %-#10.7g\n', mfilename, string(tracer_cell(idx)), ...
        tracerError(col), tracerNorm(col), xnrm)
end

end % unscrambleAndSaveParFor function
