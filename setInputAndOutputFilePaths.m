function sim = setInputAndOutputFilePaths(sim, args) % tracer_loop inputRestartFile time_step_hr recalculate_PQ_inv "short_circuit G() and phi()" logTracers
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Some tracers are known NOT work in single column with precoditioner, but
% do not crash MARBL. Use 'G()' rather than 'r(G()' as 'f()' in nsoli()

disabledPreconditoners = { 'DIC' 'ALK' 'diatC' 'spChl' 'diatChl' 'diazChl'};
disabledPreconditoners = unique( disabledPreconditoners);
sim.disabledPreconditoners = disabledPreconditoners;
fprintf('%s.m: Tracers not to be preconditioned: %s\n', mfilename, strjoin(disabledPreconditoners));


% Validate list of tracers to be solved; if one was input. If a list was
% not input, create a default one.

if  length(args) >= 1    % there is a list of tracers to be checked
    if iscell(args{1})
        sim.tracer_loop = args{1};
    else
        error('\n"tracer_loop" input arg must be a cell array even if only 1 element, {''He''}, not a %s.', class(args{1}))
    end
else                    % no list was input

    fprintf('%s.m: Creating a default ist of tracers to loop over\n', mfilename);
    tName = tracer_names(0);    % no CISO tracers
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

% sim.tracer_loop = {'spCaCO3' 'spFe' 'diazFe' };
fprintf('%s.m: *** TENATIVE *** List of tracers to loop over : %s\n', mfilename, strjoin(sim.tracer_loop));


% Tracers that we do not even try to solve.
%
% ALWAYS punt "ALT" unused methods do NOT influence other tracers, but
% waste lots of run time.

sim.excluded_tracer = { 'DIC_ALT_CO2' 'ALK_ALT_CO2' };

% These tracers do NOT work in single column; and cause "MARBL crash".
% FIXME: Maybe because preconditioner for them is messed up!

sim.excluded_tracer = [ sim.excluded_tracer ];          % FIXME???
sim.excluded_tracer = unique( sim.excluded_tracer );

% remove excluded and make sure all choices are valid...
% ...a nightmare of dealing with empty arrays and so on.
%
% boils down to finding idx into list of tracers

if numel(sim.excluded_tracer) == 0  % Nothing to exclude, move on!
    idx = [];
else
    [flag, idx] = ismember ( sim.excluded_tracer, sim.tracer_loop );
    idx = sort(idx(flag>0));
end

if any(idx)
    sim.tracer_loop([idx]) = [];
end
fprintf('%s.m: Tracers EXCLUDED from run: %s\n', mfilename, strjoin(sim.excluded_tracer));
fprintf('%s.m: *** Final *** List of tracers to loop over : %s\n', mfilename, strjoin(sim.tracer_loop));

if  length(args) >= 2    % Input restart file INCLUDING PATH
    sim.inputRestartFile = args{2};
    % FIXME: what is start_yr? it's not in saved file, but it is in fname
    sim.start_yr = 0;
else
    %     [outputArg1,outputArg2] = setInputAndOutputFilePaths(inputArg1,inputArg2);
    % sim.start_yr = 4101;  inputRestartFileStem = 'Data/InputFromAnn/restart4101.mat';
    % sim.start_yr =   0;   inputRestartFileStem = 'Data/passive_restart_init.mat'; % from netCDF 5/25/22
    sim.start_yr = 260;   inputRestartFileStem = 'Data_GP/restart_260_integrate_from_0.mat';
    % sim.start_yr = 260;   inputRestartFileStem = 'restart_260_integrate_from_0.mat';
    % sim.start_yr = 1323;  inputRestartFileStem = 'restart_0_1_output/restart_1323_DOP_sol_x1.mat';
    % sim.start_yr = 260;  inputRestartFileStem = 'restart_0_1_output/restart_260_diazFe_sol_x0.mat';
    % sim.start_yr = 260;  inputRestartFileStem = 'restart_0_1_output/restart_0_diazP_sol_x0.mat';
    % sim.start_yr = 260;   inputRestartFileStem = 'restart_0_1_output/restart_261_O2_fwd_x1.mat';

    sim.inputRestartFile = strcat(myDataDir(), inputRestartFileStem);
end
fprintf('%s.m: Input offline restart file: "%s"\n', mfilename, sim.inputRestartFile);
if ~isfile(sim.inputRestartFile)
    error("missing file or typo in name of inputRestartFile")
end
fprintf('%s.m: FIXME Output solution files will be given an arbitrary year of %d\n', mfilename, sim.start_yr);

if  length(args) >= 3    % Input time step in hours
    sim.time_step_hr = args{3};
else
    sim.time_step_hr = 3;
end

if  length(args) >= 4    % recalculate_PQ_inv?
    sim.recalculate_PQ_inv = args{4};
else
    sim.recalculate_PQ_inv = 1;
end

sim.debug_disable_phi               = 0;
sim.disable_ALL_Preconditioner      = 0;
% specify short_circuit phi()?
if  length(args) >= 5 
    sim.debug_disable_phi           = args{5};
    sim.disable_ALL_Preconditioner  = args{5};
end

if  length(args) >= 6    % log every time step?
    sim.logTracers = args{6};
else
    sim.logTracers = 1;
end

%%%%%% OUTput restart file

% Need input restart name to make our output ???
%
% "restart file" -FROM- this OFFline sim for restart on Mac, not CESM.
% Another hack is needed to move results from NK here back to CESM.

sim.yearsBetweenRestartFiles = 10;

sim.outputRestartDir = strcat(myDataDir(),'restart_0_1_output/');
fprintf('%s.m: Results will be saved in directory %s\n\n', mfilename,sim.outputRestartDir);
[status, msg, msgID] = mkdir(sim.outputRestartDir);
if status ~=1
    disp(msg); disp(msgID); disp(' ')
    keyboard
else
    clear status msg msgID
end

sim.marbl_file = 'Data/marbl_in'; % MARBL chemistry and other constants.


fprintf('%s.m: time_step_hr is %d hr\n', mfilename, sim.time_step_hr);
fprintf('%s.m: logTracers is %d\n', mfilename, sim.logTracers);
fprintf('%s.m: recalculate_PQ_inv is %d\n', mfilename, sim.recalculate_PQ_inv);
fprintf('%s.m:  disable_ALL_Preconditioner is %d\n', mfilename, sim.disable_ALL_Preconditioner);
fprintf('%s.m: debug_disable_phi is %d\n', mfilename, sim.debug_disable_phi);

% In past I debuged MARBL Carbon isotopes. "lciso_on", and that stuff, it
% probably still works but they makes everything bigger and mch slower.

sim.lciso_on = 0;   % run with Carbon Isotopes ??
sim.epsilon = sqrt(eps);
sim.logDiags = and (0, sim.logTracers) ; % Usually no diags..
sim.captureAllSelectedTracers = 0;


% % % disable all simulation, just check logic of filenames etc
% % sim.debug_disable_phi = 1;
% % sim.disable_ALL_Preconditioner = 1;
% % sim.disabledPreconditoners = []


% sim.tracer_loop = {'diatC'};
% sim.tracer_loop = fliplr(sim.tracer_loop);
% sim.tracer_loop = sim.tracer_loop(find( strcmp(sim.tracer_loop,'diatChl')):end );
sim.tracer_loop = {'DOPr' 'DONr' 'DOCr' 'O2'};
sim.tracer_loop = {'DOPr' 'DONr' 'DOCr'};
sim.tracer_loop = {'DOPr' 'DONr' 'Fe'};

tName = tracer_names(0);    % no CISO tracers

if ~all(matches(sim.tracer_loop,tName))
    errStr = sim.tracer_loop(~matches(sim.tracer_loop,tName));
    error('\n%s.m: Tracer list "sim.tracer_loop"... \n\n\t"%s"\n\n contains one or more invalid tracer names: \n\n\t"%s"', mfilename, strjoin(string(sim.tracer_loop)), strjoin(string(errStr)))
end

for idx = 1: numel(sim.tracer_loop)
    sim.tracer_loop_idx(idx) = find( strcmp(tName,sim.tracer_loop(idx)) );
end

end