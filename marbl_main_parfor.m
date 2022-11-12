function [ierr, x0_sol] = marbl_main(varargin) % tracer_loop, inputRestartFile, time_step_hr, logTracers, recalculate_PQ_inv, short_circuit
% function varargout = redplot(varargin) [varargout{1:nargout}] = plot(varargin{:},'Color',[1,0,0]); end

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

sim.forwardIntegrationOnly    = 0;      % 1 -> no NK just fwd integration
sim.num_relax_iterations      = 0;      % 0 means no relax steps, just use NK x1_sol
sim.num_forward_years         = 0;      % if fwd only, num fwd, else this inum fwd after relax step

sim.runInParallel = 1;      % parallel is hard to debug, but 2x faster
sim.verbose_debug = 1;
sim = setInputAndOutputFilePaths(sim, varargin)


% FIXME: hack in some stuff for debug
% keyboard
sim.time_step_hr = 12;

% % % disable all simulation, just check logic of filenames etc
% % sim.runInParallel = 0;
% % sim.debug_disable_phi = 1;

% % sim.disable_ALL_Preconditioner = 1;
% % sim.disabledPreconditoners = []

% % sim.recalculate_PQ_inv = 1;

% sim.tracer_loop = {'diatC'};
% sim.tracer_loop = fliplr(sim.tracer_loop);
% sim.tracer_loop = sim.tracer_loop(find( strcmp(sim.tracer_loop,'diatChl')):end );
sim
% keyboard


tName = tracer_names(0);    % no CISO tracers
% %     sim.selection = [ find( strcmp(tName,'diatC') ) ];
% %     sim.selection(ismember(sim.selection, [9,11]))=[];
% %     sim.selection = unique(sort(sim.selection));
% %     cstr = tName(sim.selection)';
% %     fprintf('%s.m: Selected tracer(s): #%d, "%s"\n', mfilename, sim.selection, string(cstr));

% if numel(sim.disabledPreconditoners)
%     ismember(tName(sim.selection), sim.disabledPreconditoners)
% end
% 
if ~all(matches(sim.tracer_loop,tName))
    errStr = sim.tracer_loop(~matches(sim.tracer_loop,tName));
    error('\n%s.m: Tracer list "sim.tracer_loop"... \n\n\t"%s"\n\n contains one or more invalid tracer names: \n\n\t"%s"', mfilename, strjoin(string(sim.tracer_loop)), strjoin(string(errStr)))
end

%%%%%%
sim.grd     = load(sim.inputRestartFile,'sim').sim.grd;
sim.domain  = load(sim.inputRestartFile,'sim').sim.domain;
if 0 && sim.debug_disable_phi
    fprintf('\n\n\t%s.m: ********* phi() is short circuited skip MTM read  *********\n\n',mfilename)
    MTM = 1;
else
    MTM = load(sim.inputRestartFile,'MTM').MTM;
end

sim = setPeek(sim);

sim.phi_years = 1;      % NK always uses 1 year integration

% for tracer_str = sim.tracer_loop
num_str = numel(sim.tracer_loop);
parfor par_idx = 1:num_str,2
    tmp_sim = sim;
    tracer_str = tmp_sim.tracer_loop (par_idx);
    parfor_inner(tmp_sim, MTM, tracer_str )
end % of loop over tracers
fprintf('...end of loop over tracers : '); toc(timer_total)


% FIXME: need to save workspace?!
logDir = strcat(sim.outputRestartDir,'/Logs/');
if ~exist(logDir, 'dir')
    mkdir(logDir)
end
save_timer = tic; fprintf('Saving (possibly) large workspace file... '); save(strcat(logDir,'last_run.mat'),'-v7.3','-nocompression'); toc(save_timer);

fprintf('... end of %s.m ', mfilename);
elapsedTime_all_loops_all_tracers = toc(timer_total);
disp(' ');
fprintf('\n%s.m: Finished outer solution loops over %d tracers\n', mfilename, numel(sim.tracer_loop));
disp(['Runtime: ', num2str(elapsedTime_all_loops_all_tracers, '%1.0f'),' (s) or ', num2str(elapsedTime_all_loops_all_tracers/60, '%1.1f'), ' (m)'])
fprintf('%s.m: Finished at %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));
end
