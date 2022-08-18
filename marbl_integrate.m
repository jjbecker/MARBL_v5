% Forward integrate MARBL 1000 years
%
%   matlab -nodisplay -nodesktop -nosplash -noFigureWindows -logfile batch.txt < marbl_integrate.m &
%
% MARBL code being used is https://github.com/marbl-ecosys/MARBL
%
% function marbl()
% "main" of cyclostationary transport version of MARBL,

fprintf('Start of %s.m: %s\n', mfilename, datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')));

% Need this to clear the "persistent" variables in "G()", "time_step()" and "calculate_forcing()"             
clear all; 

dbstop if error
format short g

addpath('MEX');
addpath(genpath('utils'));
addpath(genpath('plotting'));

% % FIXME: diary behavior is such that if renamed it's still active diary!!
diary off; diary off; diary on; diary off; diary on; diary on;

% Clean up for threads, more complex and much slower than expected.

timer_total = tic;

% killAndCleanThreads();  % Leftover threads and dumps cause BIG trouble
% toc(timer_total)

%%%%%%

% Setup the big picture parts of a simulation and/or NK solution.
% FIXME: Someday, when we know what inputs need to be, put all this a file

tName = tracer_names(0);    % no CISO tracers
selection = [ ...
    %     find(strcmp(tName,'SiO3')),...
    %     find(strcmp(tName,'DIC')),...
    find( strcmp(tName,'O2') ) ];
forwardIntegrationOnly = 1;
ck_years = 1000;
time_step_hr = 3;
logTracers = 0;

yearsBetweenRestartFiles = 10;

captureAllSelectedTracers = 0;

% DEBUG stuff
logTracers = 1;
ck_years = 3;
% time_step_hr = 12; % FAST debug
% yearsBetweenRestartFiles = 1;
% % time_step_hr = 6912/60/60;

% captureAllSelectedTracers = 1;

%%%%%%

marbl_file = 'Data/marbl_in'; % MARBL chemistry and other constants.

%%%%%% INput restart file

% start_yr = 0; inputRestartFile = 'Data/passive_restart_init.mat'; % from netCDF 5/25/22
start_yr =  70; inputRestartFile = 'Data_GP/restart_70_integrate_from_0.mat';
start_yr = 260; inputRestartFile = 'Data_GP/restart_260_integrate_from_0.mat';
% start_yr = 4101; inputRestartFile = 'Data/InputFromAnn/restart4101.mat';

fprintf('%s.m: Reading OFFline input restart file with tracers and transports: %s\n', mfilename, inputRestartFile);
load(inputRestartFile,'sim','MTM');

% We just over wrote sim struct, but now we can save the file stuff in it.
% Need input restart name to make our output...

sim.forwardIntegrationOnly  = forwardIntegrationOnly ;
sim.inputRestartFile        = inputRestartFile; 
sim.start_yr                = start_yr;
sim.selection               = selection;
sim.captureAllSelectedTracers = captureAllSelectedTracers;

sim.logTracers = logTracers;
sim.logDiags   = and (0, sim.logTracers) ; % Usually no diags..

% FIXME: lots of old code floating around...
clear inputRestartFile start_yr selection captureAllSelectedTracers logTracers 

sim.checkNeg = 0;

%%%%%% OUTput restart file

% "restart file" -FROM- this OFFline sim, if for restart on Mac, not CESM.
% Another hack is needed to move results from NK here back to CESM.

sim.outputRestartDir = myRestartDir(ck_years);
disp(['Results will be saved in directory ', sim.outputRestartDir]); disp (' ');
[status, msg, msgID] = mkdir(sim.outputRestartDir); 
if status ~=1
    disp(msg); disp(msgID); disp(' ')
    keyboard
else
    clear status msg msgID
end


sim.yearsBetweenRestartFiles = yearsBetweenRestartFiles;       % YEARS between drops of restart

%%%%%%

% parallel is hard to debug, but 2x faster

sim.runInParallel = 1;
if (sim.runInParallel)
    sim.number_of_threads = 4; end

%%%%%%

% In past I debuged MARBL Carbon isotopes. "lciso_on", and that stuff, it
% probably still works but they makes everything bigger and mch slower.

sim.lciso_on = 0;   % run with Carbon Isotopes ??

%%%%%%

sim.epsilon = -sqrt(eps);

%%%%%%

% I think it's impossible to debug anything with out looking at time step
% data. If you can debug without intermediate results just turn off
% logging.
% %
% However MARBL diagnostics, while occasionally invaluable are HUGE and
% slow the sim down by at least a factor 2x, require huge memry etc; so I
% very rarely capture them...

if (or (sim.logDiags, sim.logTracers))
    disp('Setting up to log tracers and possibly MARBL diags...')
    % need grid dimensions for many things. Starting with possible log files...
    [~, sim.domain.wet_loc, sim.domain.iwet_FP, ~, sim.domain.bottom_lvl, ~]...
        = createIwetEtc(sim.domain.M3d);

    % Possibly Record everything at a single location, water level. But which?
    % FIXME: need a clean way to convert (lat,lon) to (iLat,iLon), but in the
    % mean time; iterate thise by hand...
    %
    % iLat = 50; iLon = 81; iLvl = 10;    % Galapagos(0.3N, 108.7W)   iFp = 6496
    % iLat = 77; iLon = 60; iLvl = 10;    % NorPac  (45.7N, 176.8E)   iFp = 4629
    % iLat = 50; iLon = 61; iLvl = 1;     % Dateline  ( 0.3N,  179.3E)   iFp = 4702
    % iLat = 49; iLon = 11; iLvl = 10;    % Zulu =  ( 0.3N, 0.7E)     iFp = 1049
    % iLat = 58; iLon = 50; iLvl = 10;    % Palau = ( 5.6N, 139.7E)   iFp = 3704
    % iLat = 50; iLon = 28; iLvl = 10;    % IO      ( 0.3N,  50.5E)   iFp = 2080
%     iLat = 57; iLon =  3; iLvl = 10;    % AF 447 =  ( 4.7N, -29.5E)     iFp = 1049
%     iLat = 20; iLon =  95; iLvl = 4;    % "-48" =  ( -45.695N, -58.3E)     iFp = 31045 iCol 7462
    iLat = 22; iLon =  97; iLvl = 5;    % "+44" =  ( -40.425N, -51.1E)     iFp = 31045 icol 7587
    % Check that! Make a map!
    % first get iFp on level 1, Simpy put: on level 1, iFp = iCol...

    iFp = coordTransform_xyz2fp(iLat, iLon, 1, sim);
    [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(iFp, sim, 999); title('Time Series Localtion')
%     [~, ~, ~, ~, ~, ~] = coordTransform_fp2xyz(iFp, sim); 
    sim.time_series_loc = iFp ;
    % % % ... now set the level
    sim.time_series_lvl = iLvl;
    disp(['Time series(loc,lvl) = (', num2str(sim.time_series_loc), ', ', num2str(sim.time_series_lvl),')']);
    [~,~,~, lat, lon, ~] = col2latlon(sim, sim.time_series_loc);
    disp(['Time series(lat,lon, depth) = (', num2str(lat,'%.1f'), ' N, ', num2str(lon,'%.1f'),' E, ',num2str(sim.grd.zt(sim.time_series_lvl),'%.1f'),' m))'])
    disp(' ')
    clear iLat iLon iLvl lat lon iFp
else
    % avoid messy code in parallel;
    % just set a default legal array idx
    sim.time_series_loc = 1 ;
    sim.time_series_lvl = 1;
end

%%%%%% End of "inputs"
toc(timer_total)

% Stuff below is not a simulation input. It is code to setup grids, etc a
% forward integration, or NK(), or...
%
[sim, bgc, ~, time_series, forcing] = init_sim(marbl_file, sim.inputRestartFile, sim, ck_years, time_step_hr);
if (or (sim.logDiags, sim.logTracers))
    getMemSize(time_series,1e3);
end

% DEBUG: run on a few locations, rather global, MUCH faster debug something
% sim.domain.num_wet_loc = 1; % comment out too run entire world

% scale units to unity, using values of global mean of initial
% need a column vector

sim = calc_global_moles_and_means(bgc, sim);

%%%%%%

toc(timer_total)

% =============== This is a forward "test" integration ================

if (sim.checkNeg) % check for negatives in init file
    c0 = bgc2nsoli(sim, bgc.tracer);
    fprintf('Fraction of tracers(:) <0 = %s\n', num2str(sum(c0<0)/numel(c0)));
    find(c0<0);
    [negativesFound] = negative_tracer_catcher(sim,bgc);
    fprintf('Found %d (%.0f%%) negatives in %.2e tracers in %s\n', ...
        negativesFound, 100*negativesFound/numel(bgc.tracer(:)), numel(bgc.tracer(:)), sim.inputRestartFile)
    %     keyboard
end

if (sim.forwardIntegrationOnly) % if just running forward sim

    [sim, bgc, time_series] = phi(sim, bgc, time_series, forcing, MTM);
    %     [negativesFound] = negative_tracer_catcher(sim,bgc);
    %     fprintf('Found %d (%.0f%%) negatives in %.2e tracers after %d year integration\n', ...
    %         negativesFound, 100*negativesFound/numel(bgc.tracer(:)), numel(bgc.tracer(:)),ck_years)

    % save final workspace and diary...

    years_gone_by = round(sim.num_time_steps *sim.dt /sim.const.sec_y);

    if (1) % dumping everything at end of sim
        allFile = sprintf('%s/all_%d_done.mat', sim.outputRestartDir, round(years_gone_by));
        fprintf('%s.m: saving %s\n', 'marbl_integrate', allFile);
        save(allFile,'-v7.3');
        movefile ('diary', sim.outputRestartDir); diary off; diary on; toc(timer_total);
    end
end

