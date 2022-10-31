function [sim, bgc_struct] = init_marbl (marbl_file, sim, bgc_struct, surf_forcing)

% read settings file of choice (if any)

% fprintf('%s.m: using arbitrary default tracers...\n',mfilename);
% mfilename() is name of this file...
ioerr = mex_marbl_driver('read_settings_file', marbl_file);
if (ioerr), error('read_settings_file'), end

% FIXME: writing chemistry file fails because varcount = 0 in F90 code.
% ioerr = mex_marbl_driver('write_settings_file', 'new.input');
%
% FIXME: something as simple as changing a setting... sigh
% FIXME: Must hack out "ciso_on" line in defaults file if we set it here...
% FIXME: should read file or better ask MARBL for value of sim.lciso_on
if sim.lciso_on
    mex_marbl_driver('put setting', 'ciso_on = .true.');
end

% initialize a single MARBL instance with a single column in it so we can
% read its dimension, name, and unit

[ marbl_log, ~, ...
    ~, diag_cnt, ...
    ~ , surf_diag_cnt ] ...
    = mex_marbl_driver ( 'init', ...
    sim.domain.dzt *sim.domain.MARBL_depth_per_m, ...
    sim.domain.zw  *sim.domain.MARBL_depth_per_m, ...
    sim.domain.zt  *sim.domain.MARBL_depth_per_m );

firstErr = 0;
for i = 1:size(marbl_log,1)
    if (strfind(marbl_log(i,:),'ERROR')) >0
        if (firstErr == 0)
            firstErr = i;
        end
    end
end
if (firstErr >0)
    disp(marbl_log(firstErr,:));
    keyboard
end
mex_marbl_driver('set_depth', size(sim.domain.M3d,3) )

% initialze forcing first so we can use CISO values to init tracers
% These are make believe values at one water columns

% [bgc_struct.value.surf_forcing, bgc_struct.value.forcing] ...
%     = default_forcings ( sim, bgc_struct );

% mex_marbl_driver ( 'restore_surface_flux_forcings',      bgc.surf_forcing );
% mex_marbl_driver ( 'restore_interior_tendency_forcings', bgc.forcing );


% initialize everything in MEX to something very roughly like SMOW or it
% crashes on first calls...

[~, bgc_struct.value.state ] = init_states();
[bgc_struct.value.tracer]    = init_tracers  ( surf_forcing, sim.lciso_on );
% [bgc_struct.value.old_tracer]= init_tracers  ( surf_forcing, sim.lciso_on );

% No initialization needed for tendency, diags, and fluxes which are MARBL
% outputs but we read them to get dimensions, then pre-allocate log, so we
% can save results quickly...

bgc_struct.value.tendency  = nan(size(mex_marbl_driver ( 'interior_tendencies' )))';
bgc_struct.value.surf_flux = nan(size(mex_marbl_driver ( 'surface_fluxes' )));
% FIXME: these are hacked to be matrix, but diags in MARBL are complex structs...
bgc_struct.value.diag      = nan(size(mex_marbl_driver ( 'interior_tendency_diags')))';
bgc_struct.value.surf_diag = nan(size(mex_marbl_driver ( 'surface_flux_diags')))';

bgc_struct.name.forcing      = interior_forcing_names ();
bgc_struct.name.surf_forcing = surface_forcing_names ( sim.lciso_on );
bgc_struct.name.river_flux   = river_flux_names ( sim.lciso_on ); % FIXME: flux ~= tend, but tend is what we eventually need
bgc_struct.name.tracer       = tracer_names ( sim.lciso_on );
% bgc_struct.name.old_tracer   = tracer_names ( sim.lciso_on );
bgc_struct.name.state        = state_names  ();
bgc_struct.name.tendency     = tracer_names ( sim.lciso_on );
bgc_struct.name.sfo          = surface_flux_output_names ();
bgc_struct.name.diag         = diag_names   ( 'interior_tendency_diags', diag_cnt )';
bgc_struct.name.surf_diag    = diag_names   ( 'surface_flux_diags'     , surf_diag_cnt )';
bgc_struct.name.surf_flux    = tracer_names ( sim.lciso_on );

bgc_struct.unit.tracer       = tracer_units ( sim.lciso_on );
% bgc_struct.unit.old_tracer   = tracer_units ( sim.lciso_on );
bgc_struct.unit.forcing      = interior_forcing_names ();
bgc_struct.unit.surf_forcing = surface_forcing_names  (sim.lciso_on);
bgc_struct.unit.river_flux   = river_flux_units(sim.lciso_on);

end % init_marbl.m
