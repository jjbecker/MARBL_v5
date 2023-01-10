function sim = setInputAndOutputFilePaths(args) % tracer_loop inputRestartFile time_step_hr recalculate_PQ_inv "short_circuit G() and phi()" logTracers
% function sim = setInputAndOutputFilePaths(sim, args) % tracer_loop inputRestartFile time_step_hr recalculate_PQ_inv "short_circuit G() and phi()" logTracers
%UNTITLED Summary of this function goes here
%   args = {a, b, c}, and a, b, and c can be cells, strings, or ints
%
% tmpTracer_loop  = {'zooC' 'spC'}; disp(tmpTracer_loop)
% tmpInputFile    = strcat(myDataDir(), 'Data_GP/restart_260_integrate_from_0.mat'); disp(tmpInputFile)
% tmpTime_step_hr = 12; disp(tmpTime_step_hr)
% sim = setInputAndOutputFilePaths({tmpTracer_loop, tmpInputFile, tmpTime_step_hr});
%
%   1 -  list of cells for tracers {''He''}}
%   2 - input restart file INCLUDING PATH
%   3 - time step as integer (3 or 12 work)
%   4 - recalculate_PQ_inv (integer 0 or 1)
%   5 - debug_disable_phi (integer 0 or 1)
%   6 - log file and plots (integer 0 or 1)


tName = tracer_names(0);    % no CISO tracers

% Some tracers are known NOT work in single column with precoditioner, but
% do not crash MARBL. Use 'G()' rather than 'r(G()' as 'f()' in nsoli()

disabledPreconditoners = { 'DIC' 'ALK' 'diatC' 'spChl' 'diatChl' 'diazChl'};
% FIXME; Fe might not precondition correctly
disabledPreconditoners(end+1) = {'Fe'};

disabledPreconditoners = unique( disabledPreconditoners);
sim.disabledPreconditoners = disabledPreconditoners;
fprintf('%s.m: Tracers not to be preconditioned: %s\n', mfilename, strjoin(disabledPreconditoners));

%%%
% Validate list of tracers to be solved; if one was input. If a list was
% not input, create a default one.

if  length(args) >= 1    % List of tracers to be checked
    if iscell(args{1})
        if sum(~ismember(args{1},tName)) >0
            fprintf('\n%s.m: MARBL tracers are: %s\n', mfilename, strjoin(tName))
            fprintf('%s.m: Input tracers are: %s\n', mfilename, strjoin(args{1}))
            tmp = args{1};
            fprintf('\n%s.m: Bogus tracers are: %s\n', mfilename, strjoin(tmp(~ismember(args{1},tName))))
            error('One or more of tracers is NOT a MARBL tracers')
        end
        sim.tracer_loop = args{1};
    else
        error('\n"tracer_loop" input arg must be a cell array possibly of cell arrays even if only 1 element, {{''He''}}, not a %s.', class(args{1}))
    end

else                    % no list was input

    fprintf('%s.m: Creating a default ist of tracers to loop over\n', mfilename);
    sim.tracer_loop = tName;

    % Most of these work very well in single tracer solution...
    % FIXME: rather than get accurate solution, loop over all tracers and
    % try to reduce G to 1% of starting value, while looping over all
    % tracers, and then repeat until to get final very accurate result
    % where G is sqrt(eps).
    %
    % Always need a selected tracer! For plot time series, or solve!
    % sim.tracer_loop = {'DOPr' 'DONr' 'DOCr' 'O2' 'DON' 'DOC' 'DOP' 'diatSi' 'spCaCO3' 'diazFe' };
    % sim.tracer_loop = {'DOC' 'DOP' 'spCaCO3' 'diatSi' 'diazFe' };
    % sim.tracer_loop = {'spCaCO3' 'spFe' 'diazFe' };
    % sim.tracer_loop = {'diatFe' 'spP' 'diatP' 'diazP'};
    % sim.tracer_loop = {'zooC' 'spC' 'spCaCO3' 'diatC' 'diatFe' 'diatSi' 'diazC' 'NH4' 'Fe' 'DOP' 'diatChl' 'diazChl'};

end

%%%
% sim.tracer_loop = {'spCaCO3' 'spFe' 'diazFe' };
fprintf('%s.m: *** TENATIVE *** List of tracers to loop over : %s\n', mfilename, strjoin(sim.tracer_loop));


% Tracers that we do not even try to solve.
%
% ALWAYS punt "ALT" unused methods do NOT influence other tracers, but
% waste lots of run time.
%
% These tracers do NOT work in single column; and cause "MARBL crash".
% FIXME: Maybe because preconditioner for them is messed up!
%
% remove excluded and make sure all choices are valid...
% ...a nightmare of dealing with empty arrays and so on.
%
% boils down to finding idx into list of tracers

sim.excluded_tracer = unique({ 'DIC_ALT_CO2' 'ALK_ALT_CO2' });
if numel(sim.excluded_tracer) == 0  % Nothing to exclude, move on!
    idx = [];
else
    [flag, idx] = ismember ( sim.excluded_tracer, sim.tracer_loop );
    if (length(args) >0)
        idx = sort(idx(flag>0));
    else
        idx = 0;
    end
end
if any(idx)
    fprintf('%s.m: Tracers ALWAYS EXCLUDED from run: %s\n', mfilename, strjoin(sim.excluded_tracer));
    %     error(sprintf('\n%s.m: illegal INPUT Tracers: %s\n', mfilename, strjoin(sim.tracer_loop(idx))));
    fprintf('\n%s.m: Removing illegal INPUT Tracers: %s\n', mfilename, strjoin(sim.tracer_loop(idx)));
    sim.tracer_loop([idx]) = [];
end
fprintf('%s.m: *** Final *** List of tracers to loop over : %s\n', mfilename, strjoin(sim.tracer_loop));

%%%
if  length(args) >= 2    % Input restart file INCLUDING PATH
    if ischar(args{2})
        sim.inputRestartFile = args{2};
        % FIXME: what is start_yr? it's not in saved file, but it is in fname
        %         start_yr = 260;
        start_yr = 666;
        sim.start_yr = start_yr;
        fprintf('%s.m: FIXME Output solution files will be given an arbitrary year of %d\n', mfilename, sim.start_yr);
    else
        error("inputRestartFile must be string")
    end
else
    %     [outputArg1,outputArg2] = setInputAndOutputFilePaths(inputArg1,inputArg2);
    % sim.start_yr = 4101;  inputRestartFileStem = 'Data/InputFromAnn/restart4101.mat';
    % sim.start_yr =   0;   inputRestartFileStem = 'Data/passive_restart_init.mat'; % from netCDF 5/25/22
    % sim.start_yr = 260;   inputRestartFileStem = 'Data_GP/restart_260_integrate_from_0.mat';
    % sim.start_yr = 260;   inputRestartFileStem = 'restart_0_1_output//restart_260_integrate_from_0_single_col_DOPr_DONr_Fe.mat';
    % sim.start_yr = 260;   inputRestartFileStem = 'restart_0_1_output//restart_260_integrate_from_0_single_col_DOPr_DONr_Fe_single_col_DOPr_DONr_Fe.mat'
    % sim.start_yr = 260;   inputRestartFileStem = 'restart_0_1_output//restart_260_integrate_from_0_single_col_DOPr_DONr_Fe_single_col_DOPr_DONr_Fe_single_col_DOPr_DONr_Fe.mat';

    % sim.start_yr = 260;   inputRestartFileStem = 'restart_260_integrate_from_0.mat';
    % sim.start_yr = 1323;  inputRestartFileStem = 'restart_0_1_output/restart_1323_DOP_sol_x1.mat';
    % sim.start_yr = 260;  inputRestartFileStem = 'restart_0_1_output/restart_260_diazFe_sol_x0.mat';
    % sim.start_yr = 260;  inputRestartFileStem = 'restart_0_1_output/restart_0_diazP_sol_x0.mat';
    % sim.start_yr = 260;   inputRestartFileStem = 'restart_0_1_output/restart_261_O2_fwd_x1.mat';

    sim.inputRestartFile = strcat(myDataDir(), inputRestartFileStem);
end
fprintf('%s.m: Input offline restart file: "%s"\n', mfilename, sim.inputRestartFile);
if ~isfile(sim.inputRestartFile)
    error("missing file, or file path is not fully qualified, or typo in name of inputRestartFile")
end

if  length(args) >= 3    % Input time step in hours
    if isnumeric(args{3})
        fprintf('%s.m: time step is %g hr\n', mfilename, args{3})
        if sum(~ismember(args{3},[3 12])) >0
            error("Time step must be 3 or 12 (FIXME: or possibly other integer multiples of 3 because 'SOLAR_3hr_forcing.mat' is 3hr time step...")
        end
        sim.time_step_hr = args{3};
    else
        errMsg = sprintf('%s.m: time step input must be NUMERIC, not string or cell it is (possibly string):\n', mfilename)
        disp(args{3})
        error(errMsg)
    end
else
    sim.time_step_hr = 3;
end

if  length(args) >= 4    % recalculate_PQ_inv?
    if isnumeric(args{4})
        if sum(~ismember(args{4},[0 1])) >0
            error("recalculate_PQ_inv must be 0 or 1")
        end
    else
        error("recalculate_PQ_inv must be numeric")
    end
    sim.recalculate_PQ_inv = args{4};
else
    sim.recalculate_PQ_inv = 1;
end

sim.debug_disable_phi               = 0;
sim.disable_ALL_Preconditioner      = 0;
% specify short_circuit phi()? If so disable preconditiomer also...
if  length(args) >= 5
    sim.debug_disable_phi           = args{5};
    sim.disable_ALL_Preconditioner  = args{5};
end

if  length(args) >= 6    % log every time step?
    sim.logTracers = args{6};
else
    sim.logTracers = 0;
end

%%% OUTput restart file
% Need input restart name to make our output ???
%
% "restart file" -FROM- this OFFline sim for restart on Mac, not CESM.
% Another hack is needed to move results from NK here back to CESM.

sim.outputRestartDir = strcat(myDataDir(),'restart_0_1_output/');
fprintf('%s.m: Results will be saved in directory %s\n\n', mfilename,sim.outputRestartDir);
[status, msg, msgID] = mkdir(sim.outputRestartDir);
if status ~=1
    disp(msg); disp(msgID); disp(' ')
    keyboard
else
    clear status msg msgID
end

%%% MARBL constants 'marbl_file'
sim.marbl_file = 'Data/marbl_in'; % MARBL chemistry and other constants.


% fprintf('%s.m: time_step_hr is %d hr\n', mfilename, sim.time_step_hr);
% fprintf('%s.m: logTracers is %d\n', mfilename, sim.logTracers);
% fprintf('%s.m: recalculate_PQ_inv is %d\n', mfilename, sim.recalculate_PQ_inv);
% fprintf('%s.m:  disable_ALL_Preconditioner is %d\n', mfilename, sim.disable_ALL_Preconditioner);
% fprintf('%s.m: debug_disable_phi is %d\n', mfilename, sim.debug_disable_phi);

% In past I debuged MARBL Carbon isotopes. "lciso_on", and that stuff, it
% probably still works but they makes everything bigger and mch slower.

sim.lciso_on = 0;   % run with Carbon Isotopes ??
sim.epsilon = sqrt(eps);
sim.logDiags = and (0, sim.logTracers) ; % Usually no diags..
% sim.captureAllSelectedTracers = 0;  % log a years worth of tendency for preconditioer experiments

sim.runInParallel           = 0;    % parfor can't use spmd inside, at least I can not make that work
sim.verbose_debug           = 0;
sim.num_single_tracer_relax_iters = 0;    % 0 means no relax steps, just use NK x1_sol
sim.num_forward_iters       = 0;    % num of bgc = phi(bgc) loops, but sim.phi_years can be >1.

% NK always uses 1 year integration
sim.phi_years           = 1;    
sim.forwardIntegrationOnly  = 0;    % 1 -> no NK just fwd integration
if sim.forwardIntegrationOnly
    % FIXME: this is hacked in marbl_main_parfor to the number of interest.
    % if not doing NK, and just doing forward sims, ok then...
    keyboard
    sim.phi_years           = 1;    % VERY special case. >1 to check end of year bugs.
end

sim.grd     = load(sim.inputRestartFile,'sim').sim.grd;
sim.domain  = load(sim.inputRestartFile,'sim').sim.domain;
[sim, ~]    = calculate_depth_map_and_volumes(sim);

sim  = setPeek(sim);

%%%
for idx = 1: numel(sim.tracer_loop)
    sim.tracer_loop_idx(idx) = find( strcmp(tName,sim.tracer_loop(idx)) );
end

%%%
% These are control for nsoli(), that is called from marbl_solve
sim.maxfeval    = 11;       % assumes we input f(x0)
sim.rtol        = 1e-2;     % stop if norm(drift,2) < 10% of G(x0)

end % of function


%%%
% % % % % % % fprintf('\n %%%%%%%%%%%%%%%%   HACKERY !!!! %%%%%%%%%%%%%%%% \n\n');
% % % % % % %
% % % % % % % % FIXME: hack in some stuff for debug
% % % % % % % sim.time_step_hr = 12;
% % % % % % % % % % disable all simulation, just check logic of filenames etc
% % % % % % % % sim.recalculate_PQ_inv         = 0;
% % % % % % % % sim.debug_disable_phi          = 1;
% % % % % % % % sim.disable_ALL_Preconditioner = 1;
% % % % % % % % % sim.disabledPreconditoners = []
% % % % % % %
% % % % % % %
% % % % % % % % sim.tracer_loop = {'diatC'};
% % % % % % % % sim.tracer_loop = fliplr(sim.tracer_loop);
% % % % % % % % sim.tracer_loop = sim.tracer_loop(find( strcmp(sim.tracer_loop,'diatChl')):end );
% % % % % % % sim.tracer_loop = {'DOP' 'DON' 'DOC' 'O2'};
% % % % % % % sim.tracer_loop = {'DOP' 'DOC'};
% % % % % % % % sim.tracer_loop = {'DOPr' 'DONr' 'DOCr'};
% % % % % % % % sim.tracer_loop = {'DOPr' 'DONr' 'Fe' 'O2'};
% % % % % % % % sim.tracer_loop = {'Fe' 'O2'};
% % % % % % % % sim.recalculate_PQ_inv = 0;
% % % % % % %
% % % % % % % hackedTracerLoop = sim.tracer_loop
% % % % % % % if ~all(matches(sim.tracer_loop,tName))
% % % % % % %     errStr = sim.tracer_loop(~matches(sim.tracer_loop,tName));
% % % % % % %     error('\n%s.m: Tracer list "sim.tracer_loop"... \n\n\t"%s"\n\n contains one or more invalid tracer names: \n\n\t"%s"', mfilename, strjoin(string(sim.tracer_loop)), strjoin(string(errStr)))
% % % % % % % end

