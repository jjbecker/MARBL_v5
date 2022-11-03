function sim = setInputAndOutputFilePaths(sim, args) % tracer_loop inputRestartFile time_step_hr recalculate_PQ_inv "short_circuit G() and phi()" logTracers
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% % args = varargin{1};

% disp(['Number of provided inputs: ' num2str(length(args))])
% disp(['Number of requested outputs: ' num2str(nargout)])
% disp("Total number of input arguments: " + length(args))
% formatSpec = "Size of varargin cell array: %dx%d";
% str = compose(formatSpec,size(args));
% disp(str)
% celldisp(args)

% Return request data % for k = 1:nargout
%     varargout{k} = k;
% end

if  length(args) >= 1    % Input restart file
    if iscell(args{1})
        sim.tracer_loop = args{1};
    else
        error('\n"tracer_loop" input arg must be a cell array even if only 1 element, {''He''}, not a %s.', class(args{1}))
    end
else
    % always need a selected tracer! For plot time series, or solve
    % Most of these work very well in single tracer solution...
    % FIXME: rather than get accurate solution, loop over tracers and try to reduce G
    % to 1% of starting value, while looping over all tracers, and then
    % repeat until to get final very accurate result where G is sqrt(eps)
    % sim.tracer_loop = {'DOPr' 'DONr' 'DOCr' 'O2' 'DON' 'DOC' 'DOP' 'diatSi' 'spCaCO3' 'diazFe' };
    % sim.tracer_loop = {'DOC' 'DOP' 'spCaCO3' 'diatSi' 'diazFe' };
    sim.tracer_loop = {'spCaCO3' 'spFe' 'diazFe' };
    % sim.tracer_loop = {'O2' };
end
fprintf('%s.m: Loop over tracers : %s\n', mfilename, strjoin(sim.tracer_loop));

if  length(args) >= 2    % Input restart file INCLUDING PATH
    sim.inputRestartFile = args{2};
    % FIXME: what is start_yr? it's not in saved file, but it is in fname
    sim.start_yr = 0;
else
    %     [outputArg1,outputArg2] = setInputAndOutputFilePaths(inputArg1,inputArg2);
    % sim.start_yr = 4101;  inputRestartFileStem = 'Data/InputFromAnn/restart4101.mat';
    % sim.start_yr =   0;   inputRestartFileStem = 'Data/passive_restart_init.mat'; % from netCDF 5/25/22
    % sim.start_yr = 260;   inputRestartFileStem = 'Data_GP/restart_260_integrate_from_0.mat';
    sim.start_yr = 1323;  inputRestartFileStem = 'restart_0_1_output/restart_1323_DOP_sol_x1.mat';
    % sim.start_yr = 260;   inputRestartFileStem = 'restart_0_1_output/restart_261_O2_fwd_x1.mat';

    sim.inputRestartFile = strcat(myDataDir(), inputRestartFileStem);
end
fprintf('%s.m: FIXME Loop over tracers starts with year %d of offline restart file: "%s"\n', mfilename, sim.start_yr, sim.inputRestartFile);
if ~isfile(sim.inputRestartFile)
    error("missing file or typo in name of inputRestartFile")
end

if  length(args) >= 3    % Input time step in hours
    sim.time_step_hr = args{3};
else
    sim.time_step_hr = 3;
end
fprintf('%s.m: time_step_hr is %d hr\n', mfilename, sim.time_step_hr);

if  length(args) >= 4    % recalculate_PQ_inv?
    sim.recalculate_PQ_inv = args{4};
else
    sim.recalculate_PQ_inv = 1;
end

sim.debug_disable_phi = 0;
if  length(args) >= 5                    % specify short_circuit phi()?
    sim.debug_disable_phi   = args{5};
    sim.debug_PQ_inv        = sim.debug_disable_phi;
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


fprintf('%s.m: logTracers is %d\n', mfilename, sim.logTracers);
fprintf('%s.m: recalculate_PQ_inv is %d\n', mfilename, sim.recalculate_PQ_inv);
fprintf('%s.m: debug_PQ_inv is %d\n', mfilename, sim.debug_PQ_inv);
fprintf('%s.m: debug_disable_phi is %d\n', mfilename, sim.debug_disable_phi);

% In past I debuged MARBL Carbon isotopes. "lciso_on", and that stuff, it
% probably still works but they makes everything bigger and mch slower.

sim.lciso_on = 0;   % run with Carbon Isotopes ??
sim.epsilon = -sqrt(eps);
sim.logDiags = and (0, sim.logTracers) ; % Usually no diags..
sim.captureAllSelectedTracers = 0;

end